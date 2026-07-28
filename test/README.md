# TeraLCP threshold tests

Self-contained checks for TeraLCP's pfp-thresholds-style threshold output
(`-othresholds`). They need only tools this repository builds — ropebwt3 and
`src/TeraLCP/TeraLCP` — and no external dependencies at run time.

Run from the repo root:

```
make -C src/TeraLCP        # builds TeraLCP (and ropebwt3)
bash test/run_threshold_tests.sh
# or: make -C src/TeraLCP test-thresholds
```

## What is checked

For each case the runner builds the FMD with `ropebwt3 build -R -d` and then:

1. **Correctness vs an independent reference.** TeraLCP's `.thr`/`.thr_pos` must
   match the committed goldens, which were produced **by pfp-thresholds** (the
   reference MONI-family implementation), not by TeraLCP. Because `ropebwt3 build
   -R -d` produces the same BWT pfp built for these inputs, an exact byte match is
   expected.
2. **Path agreement.** The fast parallel path (default) and the non-destructive
   phi-walk path (used when `-orlcp` is also requested) must produce identical
   `.thr`/`.thr_pos`.
3. **`-f rlbwt` round-trip.** Feeding TeraLCP's own `.bwt.heads`/`.bwt.len` back
   through the generic run-length-BWT input must reproduce the `.thr`/`.thr_pos`.

## Cases / fixtures

Each `<case>/` directory holds the input FASTA and the pfp-thresholds goldens:

| file | meaning |
|------|---------|
| `<name>.fa` | input sequence(s) |
| `<name>.fa.thr` | golden threshold values (5-byte LE per BWT run), from pfp-thresholds |
| `<name>.fa.thr_pos` | golden threshold positions (5-byte LE per BWT run), from pfp-thresholds |

- `single_string` — a single short sequence.
- `gattacat` — small repetitive sequence (117 runs).
- `shred1_mini` — shredded reads (~14.9k runs), a larger non-trivial case.

## Broader validation

The grlBWT-to-TeraLCP adapter and the end-to-end cross-tool validation against Movi
— a Movi index built from TeraLCP's thresholds answers queries identically to the
pfp-thresholds pipeline, over both the ropeBWT3 and grlBWT construction paths — live
in [thresh-tools](https://github.com/BenLangmead/thresh-tools).
