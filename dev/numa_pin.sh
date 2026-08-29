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
# both arms always share the policy. Whether binding is as *reproducible* as
# interleaving is not yet measured, which is the only reason the default here
# has not changed. See docs/tuning-for-your-data.md.
#
# NOTE this is a *measurement* tool. Interleaving is not a speed recommendation:
# at 64 threads it matches default placement on mean run_dada (372s vs 372s).
# What it buys is predictability.
#
# Usage:
#   source "$(dirname "$0")/numa_pin.sh"
#   numa_pin_init "$THREADS"        # sets NUMA_PREFIX (may be empty)
#   $NUMA_PREFIX some_command ...
#
# Env overrides:
#   NUMA_POLICY=interleave  (default) numactl --interleave=all
#   NUMA_POLICY=bind        numactl --cpunodebind=$NUMA_NODE --membind=$NUMA_NODE
#                           — lowest-latency for serial scattered phases, but
#                             slows the parallel map; refuses to oversubscribe
#   NUMA_POLICY=none        no numactl — use when comparing against history
#                             measured under default placement
#   NUMA_NODE=0             which domain, for NUMA_POLICY=bind
#   NUMA_PIN=0              legacy alias for NUMA_POLICY=none
# ---------------------------------------------------------------------------

# Physical cores in NUMA node $1 — distinct core ids, so SMT siblings are not
# double counted. Empty if the topology cannot be read.
_numa_cores_in_node() {
    command -v lscpu >/dev/null 2>&1 || return 0
    lscpu -p=CPU,CORE,NODE 2>/dev/null \
        | grep -v '^#' \
        | awk -F, -v n="$1" '$3==n {print $2}' \
        | sort -u | wc -l | tr -d ' '
}

# Sets NUMA_PREFIX (possibly empty) and explains the decision on stderr.
# $1 = thread count the benchmark will use.
numa_pin_init() {
    local threads="${1:-0}"
    local policy="${NUMA_POLICY:-interleave}"
    local node="${NUMA_NODE:-0}"
    NUMA_PREFIX=""

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

    local candidate
    case "$policy" in
        interleave)
            candidate="numactl --interleave=all"
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
            candidate="numactl --cpunodebind=$node --membind=$node"
            ;;
        *)
            echo "[numa] unknown NUMA_POLICY=$policy (expected interleave|bind|none) — default placement" >&2
            return 0
            ;;
    esac

    # numactl can exist but be forbidden (containers, restricted cpusets).
    if ! $candidate true >/dev/null 2>&1; then
        echo "[numa] numactl present but '$candidate' failed (container or restricted cpuset?) — default placement" >&2
        return 0
    fi

    NUMA_PREFIX="$candidate"
    echo "[numa] fixed placement: $NUMA_PREFIX" >&2
    echo "[numa]   for measurement reproducibility — see the header of dev/numa_pin.sh" >&2
}
