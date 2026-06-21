/* Oracle for WFA2-lib ends-free + nonzero-match-score divergence (issue #102).
 *
 * Links the *C++* WFA2-lib static lib (built by run.sh) and asks, on the
 * committed dada2-rs reproducers, whether the C++ library returns the OPTIMAL
 * ends-free alignment under DADA2 scoring (match +5) — i.e. whether #102 is
 * fixed upstream. Nothing here ships; this is a throwaway oracle (see README).
 *
 * Penalty config mirrors the pure-Rust adapter exactly:
 *   match=-5 (reward), mismatch=4, gap_opening=0, gap_extension=8, ends-free
 *   on both ends of both sequences. With match=-5 the WFA penalty score equals
 *   the DADA2 score, so a WFA result below the known optimum is a genuine
 *   optimality miss, not a penalty-space tie.
 */
#include <stdio.h>
#include <string.h>
#include "wavefront/wavefront_align.h"

/* DADA2 score of a WFA CIGAR op string; leading/trailing indels are free. */
static int dada2_score(const char *ops, int n) {
    int lead = 0, trail = 0;
    while (lead < n && (ops[lead] == 'I' || ops[lead] == 'D')) lead++;
    while (trail < n && (ops[n - 1 - trail] == 'I' || ops[n - 1 - trail] == 'D')) trail++;
    int s = 0;
    for (int i = 0; i < n; i++) {
        char o = ops[i];
        if (o == 'M') s += 5;
        else if (o == 'X') s += -4;
        else if (o == 'I' || o == 'D') s += (i < lead || i >= n - trail) ? 0 : -8;
    }
    return s;
}

static int run(const char *name, const char *p, const char *t, int optimal) {
    int pl = strlen(p), tl = strlen(t);
    wavefront_aligner_attr_t a = wavefront_aligner_attr_default;
    a.distance_metric = gap_affine;
    a.affine_penalties.match = -5;       /* reward — the #102 trigger */
    a.affine_penalties.mismatch = 4;
    a.affine_penalties.gap_opening = 0;  /* linear gap of 8 == DADA2 gap_p */
    a.affine_penalties.gap_extension = 8;
    a.alignment_form.span = alignment_endsfree;
    a.alignment_form.pattern_begin_free = pl;
    a.alignment_form.pattern_end_free = pl;
    a.alignment_form.text_begin_free = tl;
    a.alignment_form.text_end_free = tl;
    wavefront_aligner_t *wf = wavefront_aligner_new(&a);
    wavefront_align(wf, p, pl, t, tl);
    cigar_t *c = wf->cigar;
    int n = c->end_offset - c->begin_offset;
    const char *ops = c->operations + c->begin_offset;
    int sc = dada2_score(ops, n);
    int pass = (sc == optimal);
    printf("[%-24s] opt=%d  C++=%d  %s\n    cigar=%.*s\n",
           name, optimal, sc, pass ? "PASS (optimal)" : "FAIL (suboptimal)", n, ops);
    wavefront_aligner_delete(wf);
    return pass;
}

int main(void) {
    int ok = 0, tot = 0;
    /* Mode 1 — free leading/trailing end-gap crediting (the headline #102 symptom). */
    tot++; ok += run("Mode1 leading (255)",
        "AACAGCGCAAACCAACTCGCTAGCTAGCAAAATCTTGTGTTTCTGCCTAGCG",
        "ACAGCGCAAACCAACTCGCTAGCTAGCAAAATCTTGTGTTTCTGCCTAGCG", 255);
    tot++; ok += run("#102 thread GCGG/CG (10)", "GCGG", "CG", 10);
    tot++; ok += run("Mode1 leading (275)",
        "AACTAGACATCAAAGCCACAGAATTTGGCAGGAGATACAAGAAGTACTCCAGCTC",
        "AAACTAGACATCAAAGCCACAGAATTTGGCAGGAGATACAAGAAGTACTCCAGCTC", 275);
    /* Mode 2 — mismatch + free-end indel preferred over higher-scoring interior gap. */
    tot++; ok += run("Mode2 trailing (272)",
        "TACACGCTACCTTGGAGACCTCCTGTGCTTGCAGCCCATCAATCTTTTTACACAAGG",
        "TACACGCTACCTTGGAGACCTCCTGTGCTTGCAGCCCATCAATCTTTTTACACAGG", 272);
    tot++; ok += run("Mode2 near-end (347)",
        "TAACGACGAACGAGTTTCTAAACCAAAGAACGTAACAAGGAGTTGCCTCTTAGCCTTGGCCTACGCAGGAA",
        "TAACGACGAACGAGTTTCTAAACCAAAGAACGTAACAAGGAGTTGCCTCTTAGCCTTGGCCTACGCAGGGAA", 347);
    printf("\n%d/%d optimal.  PASS = #102 fixed for that case; FAIL = still suboptimal upstream.\n", ok, tot);
    return ok == tot ? 0 : 1;
}
