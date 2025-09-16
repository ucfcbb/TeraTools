#!/usr/bin/env bash

set -x
set -e

if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <input_file>"
    exit 1
fi

INPUT_FILE="$1"


ROPEBWT3="../../ropebwt3/ropebwt3"
BUILDER="../builder/builder"
REPEATS="./repeats"
BRUTE_FORCE_CPP="../../less_efficient_but_accurate_algos/SMEM/Repeat/brute_force_repeats.cpp"
BRUTE_FORCE_BIN="../../less_efficient_but_accurate_algos/SMEM/Repeat/brute_force_repeats"

ROPEBWT3_OUT="../../data/test.fmd"
REPEATS_OUT="rp1.out"
BRUTE_OUT="rp2.out"


echo "Building ropebwt3 index..."
$ROPEBWT3 build -L -R "$INPUT_FILE" -d -o "$ROPEBWT3_OUT"


if [[ ! -x "$BRUTE_FORCE_BIN" ]]; then
    echo "Compiling brute_force_repeats.cpp..."
    g++ -O2 -std=c++11 "$BRUTE_FORCE_CPP" -o "$BRUTE_FORCE_BIN"
fi

echo "Running builder (optbwtrl) on $ROPEBWT3_OUT..."
$BUILDER "$ROPEBWT3_OUT" "$ROPEBWT3_OUT"

echo "Running repeats.cpp on $$ROPEBWT3_OUT.optbwtrl..."
$REPEATS LM "$ROPEBWT3_OUT.optbwtrl" > "$REPEATS_OUT"

echo "Running brute_force_repeats on $ROPEBWT3_OUT..."
$BRUTE_FORCE_BIN < "$INPUT_FILE" > "$BRUTE_OUT"


tmpfile1=$(mktemp)
tail -n +5 "$REPEATS_OUT" | sed '$d' > "$tmpfile1"
header1=$(head -n 1 "$tmpfile1")
header2=$(head -n 1 "$BRUTE_OUT")

echo header1: "$header1"
echo header2: "$header2"
if [[ "$header1" != "$header2" ]]; then
    echo "Headers do not match."
    rm "$tmpfile1"
    exit 1
fi

rm -f file1.sorted file2.sorted

tail -n +2 "$tmpfile1" | sort -n -k1,1 -k2,2 -k3,3 -k4,4 > file1.sorted
tail -n +2 "$BRUTE_OUT" | sort -n -k1,1 -k2,2 -k3,3 -k4,4 > file2.sorted

echo "$header1"
diff -u file1.sorted file2.sorted

rm file1.sorted file2.sorted "$tmpfile1" "$OPT_OUT" "$REPEATS_OUT" "$BRUTE_OUT" "$ROPEBWT3_OUT"