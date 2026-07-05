#!/bin/bash
# Self-contained tests for TeraLCP's Movi-compatible threshold output.
#
# Requires only tools this repository builds: ropebwt3 (cloned/built by the
# TeraLCP Makefile) and src/TeraLCP/TeraLCP. No external tools are needed at run
# time. From the repo root:  make -C src/TeraLCP && bash test/run_threshold_tests.sh
#
# For each tiny case it checks three things:
#   (1) Correctness vs an independent reference: TeraLCP's .thr/.thr_pos must match
#       committed goldens produced by pfp-thresholds (see README.md). ropebwt3
#       'build -R -d' yields the same BWT pfp built, so an exact byte match holds.
#   (2) Path agreement: the fast parallel path (default) and the non-destructive
#       phi-walk path (taken when -orlcp is also requested) produce identical output.
#   (3) -f rlbwt round-trip: feeding TeraLCP's own .bwt.heads/.bwt.len back through
#       the generic run-length-BWT input reproduces the .thr/.thr_pos.
#
# Tiny inputs are run single-threaded (OMP_NUM_THREADS=1): TeraLCP currently has a
# known issue when the OpenMP thread count exceeds the number of BWT runs.

set -u
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RB3="${REPO_ROOT}/ropebwt3/ropebwt3"
TERALCP="${REPO_ROOT}/src/TeraLCP/TeraLCP"
export OMP_NUM_THREADS=1

CASES=(
    "single_string:single_string.fa"
    "gattacat:gattacat.fa"
    "shred1_mini:minishred1_20_002.fa"
)

for bin in "$RB3" "$TERALCP"; do
    if [[ ! -x "$bin" ]]; then
        echo "ERROR: $bin not found; build it first (e.g. make -C src/TeraLCP)." >&2
        exit 2
    fi
done

failures=0
for spec in "${CASES[@]}"; do
    dir="${spec%%:*}"; fa="${spec##*:}"
    fa_path="${REPO_ROOT}/test/${dir}/${fa}"
    gthr="${fa_path}.thr"; gpos="${fa_path}.thr_pos"
    tmp="$(mktemp -d)"
    ok=1

    "$RB3" build -R -d -o "${tmp}/index.fmd" "$fa_path" 2>/dev/null \
        || { echo "FAIL ${dir}: ropebwt3 build failed"; failures=$((failures+1)); rm -rf "$tmp"; continue; }

    # (1) pfp-thresholds-style 5-byte output (parallel path) vs pfp goldens. The default
    #     -othresholds format is the self-describing TLTHR v1; --thr-pfp emits the
    #     pfp-thresholds-style 5-byte records the committed goldens use.
    "$TERALCP" -f fmd -i "${tmp}/index.fmd" -t "${tmp}/t" -othresholds "${tmp}/A" --thr-pfp -v quiet 2>/dev/null
    cmp -s "${tmp}/A.thr"     "$gthr" || { echo "FAIL ${dir}: .thr differs from golden";     ok=0; }
    cmp -s "${tmp}/A.thr_pos" "$gpos" || { echo "FAIL ${dir}: .thr_pos differs from golden"; ok=0; }

    # (2) phi-walk path (forced by also requesting -orlcp) must equal the parallel path
    "$TERALCP" -f fmd -i "${tmp}/index.fmd" -t "${tmp}/t" -othresholds "${tmp}/B" --thr-pfp -orlcp "${tmp}/junk" -v quiet 2>/dev/null
    cmp -s "${tmp}/A.thr"     "${tmp}/B.thr"     || { echo "FAIL ${dir}: parallel vs phi-walk .thr";     ok=0; }
    cmp -s "${tmp}/A.thr_pos" "${tmp}/B.thr_pos" || { echo "FAIL ${dir}: parallel vs phi-walk .thr_pos"; ok=0; }

    # (3) -f rlbwt round-trip: TeraLCP's own heads/len reproduce the thresholds
    "$TERALCP" -f rlbwt -i "${tmp}/A" -t "${tmp}/t" -othresholds "${tmp}/C" --thr-pfp -v quiet 2>/dev/null
    cmp -s "${tmp}/A.thr"     "${tmp}/C.thr"     || { echo "FAIL ${dir}: -f rlbwt round-trip .thr";     ok=0; }
    cmp -s "${tmp}/A.thr_pos" "${tmp}/C.thr_pos" || { echo "FAIL ${dir}: -f rlbwt round-trip .thr_pos"; ok=0; }

    if [[ $ok -eq 1 ]]; then echo "PASS ${dir}"; else failures=$((failures+1)); fi
    rm -rf "$tmp"
done

if [[ $failures -gt 0 ]]; then echo "${failures} case(s) FAILED"; exit 1; fi
echo "All threshold tests passed."
