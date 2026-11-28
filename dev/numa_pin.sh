#!/usr/bin/env bash
# numa_pin.sh — sourceable helper that sets $NUMA_PREFIX for benchmark runs.
# ---------------------------------------------------------------------------
# WHY: on a multi-NUMA node, a default run lets the kernel place pages by first
# touch, and that placement re-rolls every run. On the HPCBio 128-core node
# (2 domains x 64 cores, remote distance 32 vs local 10) that made replicates of
# the *same binary* disagree by far more than most changes worth testing.
#
# Reproducibility between two replicates of one arm, 64 threads, PacBio:
#
#            run_dada   b_shuffle   store
#   default     5.8%      12.9%     22.9%
#   bound       0.9%       5.6%      1.6%
#   interleave  0.13%      0.31%     0.0%     <- default policy here
#
# That is the difference between being able to measure a 5% change and not. The
# rolling-d16 DP change (#127) was twice declared "flat, below the noise floor"
# from default-placement runs, then measured cleanly at -3.4 to -6.2% under
# fixed placement, replicates agreeing to 0.4 percentage points.
#
# WHY INTERLEAVE RATHER THAN BIND: binding to one domain is the most obvious
# fix and the wrong one. It forces every thread to draw bandwidth from a single
# node's memory controllers, which makes the *parallel* phase 21% slower
# (map 286.8s bound vs 235.4s interleaved at 64 threads; DP kernel 149,747 vs
# 119,776 ns/comp). That distorts the phase ratios these benchmarks exist to
# report. Interleaving lands within 3% of default placement while being the most
# reproducible of the three, so it measures something close to production. It
# also sets memory policy only, never CPU affinity, so unlike binding it cannot
# silently oversubscribe a thread count sized for the whole node.
#
# THREAD-COUNT CONDITIONAL -- READ THIS BEFORE APPLYING THE ABOVE (#152).
# The 21% figure was measured at 64 threads, where 64 threads bound to a 64-core
# domain saturates that domain's controllers. At 48 threads, which leaves
# headroom, binding is the FASTER policy by a wide margin -- and on both pools
# tested, not just one:
#
#   run_dada, 48 threads, interleave -> cpunodebind        ITS2 R1  ITS2 R2   16S R1
#     one job, --interleave=all                             197.8s   259.7s   708.1s
#     one job, --cpunodebind=0 --membind=0                  143.9s   185.2s   528.1s
#                                                            -27%     -29%     -25%
#
# Every phase gains, including the single-threaded ones (16S: map -27%, store
# -22%, shuffle -21%, p-update -31%), which is memory *latency* rather than
# bandwidth. Outputs byte-identical across all arms.
#
# So: do not read "interleave beats bind" as general. It holds when the thread
# count fills or exceeds a domain; it reverses when the count fits inside one.
# Benchmark numbers gathered under this helper are therefore CONSERVATIVE in
# absolute terms at sub-domain thread counts -- A/B deltas are unaffected, since
# both arms always share the policy. See docs/tuning-for-your-data.md.
#
# DEFAULT IS NOW auto (#362): bind when the thread count fits in one node's
# allowed physical cores with a quarter of the node to spare, interleave
# otherwise (48 of 64 cores binds; 64 of 64, the losing case above, does not). #360 measured the size and the
# reproducibility that kept the default at interleave. NovaSeq dada-pooled, 48
# threads, compute-5-6, job 33963, main loop (2 reps each):
#
#   --interleave=all                  983 / 968 s   (reps 1.5% apart)
#   --cpunodebind=0 --membind=0       654 / 673 s   (reps 2.9% apart)  -32%
#
# Absolute timings taken before the switch are interleave timings; compare
# against them with NUMA_POLICY=interleave.
#
# NOTE this is a *measurement* tool. Interleaving is not a speed recommendation:
# at 64 threads it matches default placement on mean run_dada (372s vs 372s).
# What it buys is predictability.
#
# Usage:
#   source "$(dirname "$0")/numa_pin.sh"
#   numa_pin_init "$THREADS"        # sets NUMA_PREFIX and NUMA_PREFIX_ARGS
#   "${NUMA_PREFIX_ARGS[@]}" some_command ...   # works under any IFS
#   $NUMA_PREFIX some_command ...               # needs the default IFS
#
# The unquoted string form relies on word splitting, so under the strict-mode
# header `IFS=$'\n\t'` it runs a command literally named "numactl --..." and
# fails; numa_pin_init warns when IFS would do that. Prefer the array.
#
# Env overrides:
#   NUMA_POLICY=auto        (default) bind if the threads fit one node's allowed
#                             physical cores and use at most 3/4 of the node,
#                             else interleave. The node is the
#                             one with the most allowed cores, then the most
#                             free memory; NUMA_NODE overrides if it fits.
#   NUMA_POLICY=interleave  numactl --interleave=all — the default until #362,
#                             for comparing against timings taken under it
#   NUMA_POLICY=bind        numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE
#                           — refuses to oversubscribe the node's cores
#   NUMA_POLICY=none        no numactl — use when comparing against history
#                             measured under default placement
#   NUMA_NODE=0             which domain, for bind (and a preference for auto)
#   NUMA_MEM=strict|preferred  memory side of a bind: --membind (strict, the
#                             default) or --preferred, which spills to another
#                             node when the bound one is full instead of
#                             failing. Strict was never slower than preferred
#                             (NovaSeq, jobs 33967/33973). Its occasional
#                             ~100 s faster screen pass is #360's two-speed
#                             screen, which strict does not cause.
#   NUMA_PIN=0              legacy alias for NUMA_POLICY=none
#
# After numa_pin_init, NUMA_RESOLVED names what was chosen: interleave,
# bind:<node>, or none.
# ---------------------------------------------------------------------------

# CPUs this process may run on (Cpus_allowed_list), one per line. Empty
# without /proc, in which case every CPU counts.
_numa_allowed_cpus() {
    [ -r /proc/self/status ] || return 0
    awk '/^Cpus_allowed_list:/ {print $2}' /proc/self/status \
        | tr , '\n' \
        | awk -F- 'NF {hi = ($2 == "") ? $1 : $2; for (i = $1; i <= hi; i++) print i}'
}

# Physical cores in NUMA node $1 that this process may use — distinct core
# ids, so SMT siblings are not double counted, and only CPUs in the cpuset, so
# a Slurm job that owns part of a socket is not credited with all of it. With
# $2 = all, every core in the node, cpuset or not. Empty if the topology
# cannot be read.
_numa_cores_in_node() {
    command -v lscpu >/dev/null 2>&1 || return 0
    local allowed=""
    [ "${2:-}" = all ] || allowed="$(_numa_allowed_cpus | paste -sd, -)"
    lscpu -p=CPU,CORE,NODE 2>/dev/null \
        | grep -v '^#' \
        | awk -F, -v n="$1" -v allowed="$allowed" '
            BEGIN { k = split(allowed, a, ","); for (i = 1; i <= k; i++) ok[a[i]] = 1 }
            $3 == n && (k == 0 || ($1 in ok)) {print $2}' \
        | sort -u | wc -l | tr -d ' '
}

# Whether auto may bind $1 threads to node $2: they must fit in the node's
# allowed cores AND leave a quarter of the whole node's cores free. Binding wins
# with headroom and loses at a full domain, where one node's memory controllers
# feed every thread: 48 of 64 cores bound was 25-32% faster than interleave
# (#152, #360), 64 of 64 was 5% slower (#127). Between the two is unmeasured.
_numa_auto_fits() {
    local threads="$1" node="$2" allowed total
    allowed="$(_numa_cores_in_node "$node")"
    total="$(_numa_cores_in_node "$node" all)"
    [ -n "$allowed" ] && [ -n "$total" ] || return 1
    [ "$threads" -le "$allowed" ] && [ $((threads * 4)) -le $((total * 3)) ]
}

# The node auto binds to for $1 threads, or nothing if none qualifies:
# NUMA_NODE if set and it qualifies, else the node with the most allowed cores,
# then the most free memory.
_numa_pick_node() {
    local threads="$1" hw n cores free best="" best_cores=0 best_free=-1
    hw="$(numactl --hardware 2>/dev/null)"
    if [ -n "${NUMA_NODE:-}" ] && _numa_auto_fits "$threads" "$NUMA_NODE"; then
        echo "$NUMA_NODE"
        return 0
    fi
    for n in $(echo "$hw" | awk '$1 == "node" && $3 == "size:" {print $2}'); do
        _numa_auto_fits "$threads" "$n" || continue
        cores="$(_numa_cores_in_node "$n")"
        free="$(echo "$hw" | awk -v n="$n" '$1 == "node" && $2 == n && $3 == "free:" {print $4}')"
        if [ "$cores" -gt "$best_cores" ] || { [ "$cores" -eq "$best_cores" ] && [ "${free:-0}" -gt "$best_free" ]; }; then
            best="$n" best_cores="$cores" best_free="${free:-0}"
        fi
    done
    echo "$best"
}

# Sets NUMA_PREFIX (possibly empty) and explains the decision on stderr.
# $1 = thread count the benchmark will use.
numa_pin_init() {
    local threads="${1:-0}"
    local policy="${NUMA_POLICY:-auto}"
    local node="${NUMA_NODE:-0}"
    local mem="${NUMA_MEM:-}"
    NUMA_PREFIX=""
    NUMA_PREFIX_ARGS=()
    NUMA_RESOLVED="none"

    # Legacy alias, kept so existing invocations keep working.
    [ "${NUMA_PIN:-}" = "0" ] && policy="none"

    if [ "$policy" = "none" ]; then
        echo "[numa] fixed placement disabled — using default first-touch" >&2
        return 0
    fi

    if ! command -v numactl >/dev/null 2>&1; then
        echo "[numa] numactl not found — default placement (expect higher run-to-run variance)" >&2
        return 0
    fi

    # Single-domain machines (and most laptops) have nothing to place.
    local nodes
    nodes="$(numactl --hardware 2>/dev/null | awk '/^available:/ {print $2}')"
    if [ -z "$nodes" ] || [ "$nodes" -le 1 ] 2>/dev/null; then
        echo "[numa] single NUMA domain — placement policy unnecessary" >&2
        return 0
    fi

    if [ "$policy" = "auto" ]; then
        local picked=""
        [ "$threads" -gt 0 ] 2>/dev/null && picked="$(_numa_pick_node "$threads")"
        if [ -n "$picked" ]; then
            echo "[numa] auto: $threads threads fit in node $picked with headroom — binding" >&2
            policy="bind" node="$picked"
        else
            if [ "$threads" -gt 0 ] 2>/dev/null; then
                echo "[numa] auto: $threads threads do not fit in 3/4 of any one node — interleaving" >&2
            else
                echo "[numa] auto: thread count not given — interleaving" >&2
            fi
            policy="interleave"
        fi
    fi

    local -a candidate
    case "$policy" in
        interleave)
            candidate=(numactl --interleave=all)
            ;;
        bind)
            # Binding restricts CPUs as well as memory, so it can turn a thread
            # count sized for the whole node into SMT contention that reads as a
            # regression in whatever is under test. Refuse rather than mislead.
            local cores
            cores="$(_numa_cores_in_node "$node")"
            if [ -n "$cores" ] && [ "$cores" -gt 0 ] 2>/dev/null && [ "$threads" -gt "$cores" ] 2>/dev/null; then
                echo "[numa] NOT binding: --threads $threads exceeds the $cores physical cores in node $node;" >&2
                echo "[numa]   binding would oversubscribe SMT siblings and corrupt the comparison." >&2
                echo "[numa]   Use NUMA_POLICY=interleave (no core restriction), or --threads <= $cores." >&2
                return 0
            fi
            case "${mem:-strict}" in
                strict) candidate=(numactl "--cpunodebind=$node" "--membind=$node") ;;
                preferred) candidate=(numactl "--cpunodebind=$node" "--preferred=$node") ;;
                *)
                    echo "[numa] unknown NUMA_MEM=$mem (expected strict|preferred) — default placement" >&2
                    return 0
                    ;;
            esac
            ;;
        *)
            echo "[numa] unknown NUMA_POLICY=$policy (expected auto|interleave|bind|none) — default placement" >&2
            return 0
            ;;
    esac

    # numactl can exist but be forbidden (containers, restricted cpusets).
    # Run as an array so a non-default IFS cannot break the probe itself.
    # Joined with explicit spaces: "${candidate[*]}" would use IFS's first char.
    local joined err
    joined="$(printf '%s ' "${candidate[@]}")"
    joined="${joined% }"
    if ! err="$("${candidate[@]}" true 2>&1)"; then
        echo "[numa] '$joined' failed — default placement. numactl said:" >&2
        echo "[numa]   ${err:-<no output>}" >&2
        return 0
    fi

    NUMA_PREFIX_ARGS=("${candidate[@]}")
    NUMA_PREFIX="$joined"
    if [ "$policy" = "bind" ]; then NUMA_RESOLVED="bind:$node"; else NUMA_RESOLVED="$policy"; fi
    echo "[numa] fixed placement: $NUMA_PREFIX" >&2
    if [ "$IFS" != $' \t\n' ]; then
        echo "[numa]   WARNING: IFS is not the default, so unquoted \$NUMA_PREFIX will not split;" >&2
        echo "[numa]   use \"\${NUMA_PREFIX_ARGS[@]}\" or restore IFS=\$' \\t\\n' around the call." >&2
    fi
    echo "[numa]   for measurement reproducibility — see the header of dev/numa_pin.sh" >&2
}
