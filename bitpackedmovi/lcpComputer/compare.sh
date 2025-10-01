#!/usr/bin/env bash


# set -x
set -e

if [[ $# -lt 4 ]]; then
    echo "Usage: $0 <num_strings_low> <num_strings_high> <len_low> <len_high> [seed]"
    exit 1
fi

NUM_STRINGS_LOW="$1"
NUM_STRINGS_HIGH="$2"
LEN_LOW="$3"
LEN_HIGH="$4"
# LEN_THRESHOLD="$6"
SEED="$5"

# Validate mode parameter
#if [[ "$MODE" != "SM" && "$MODE" != "LM" ]]; then
#    echo "Error: Mode must be either 'SM' or 'LM'"
#    exit 1
#fi

SYNTHGEN="../../data/generation/syntheticgen.cpp"
SYNTHGEN_BIN="../../data/generation/syntheticgen"
GEN_INPUT="synthetic_input.txt"

echo "Compiling syntheticgen.cpp..."
g++ -O2 -std=c++11 "$SYNTHGEN" -o "$SYNTHGEN_BIN"


if [[ -z "$SEED" ]]; then
    "$SYNTHGEN_BIN" "$NUM_STRINGS_LOW" "$NUM_STRINGS_HIGH" "$LEN_LOW" "$LEN_HIGH" > "$GEN_INPUT"
else
    "$SYNTHGEN_BIN" "$NUM_STRINGS_LOW" "$NUM_STRINGS_HIGH" "$LEN_LOW" "$LEN_HIGH" "$SEED" > "$GEN_INPUT"
fi

tmpfile="$GEN_INPUT"

INPUT_FILE=$(mktemp)
head -n -1 "$tmpfile" > "$INPUT_FILE"


ROPEBWT3="../../ropebwt3/ropebwt3"
LCPCOMPUTER="../lcpComputer/lcpComputer"
KACTL="../../less_efficient_but_accurate_algos/kactl-sa-implementation/kactl_impl"
REPEATS="./repeats"
LOADER="../loader/loader"
BRUTE_FORCE_CPP="../../less_efficient_but_accurate_algos/Repeats/brute_force_repeats.cpp"
BRUTE_FORCE_BIN="../../less_efficient_but_accurate_algos/Repeats/brute_force_repeats"

ROPEBWT3_OUT="../../data/repeatTest.fmd"
REPEATS_OUT="rp1.out"
BRUTE_OUT="rp2.out"


echo "Building ropebwt3 index..."
$ROPEBWT3 build -L -R "$INPUT_FILE" -d -o "$ROPEBWT3_OUT"


echo "Running lcpComputer on $ROPEBWT3_OUT..."
$LCPCOMPUTER "$ROPEBWT3_OUT" "$ROPEBWT3_OUT"

#echo "Validating bwt"
#$KACTL < $INPUT_FILE > bwt.out
#$LOADER "$ROPEBWT3_OUT.optbwtrl" > rbwt.out
#diff bwt.out rbwt.out
#echo "Validation successful!"
#
#
#echo "Running repeats.cpp on $ROPEBWT3_OUT.optbwtrl..."
#if [[ -n "$LEN_THRESHOLD" ]]; then
#    $REPEATS "$MODE" "$ROPEBWT3_OUT.optbwtrl" "$LEN_THRESHOLD" > "$REPEATS_OUT"
#else
#    $REPEATS "$MODE" "$ROPEBWT3_OUT.optbwtrl" > "$REPEATS_OUT"
#fi
#
#echo "Running brute_force_repeats on $INPUT_FILE..."
#if [[ -n "$LEN_THRESHOLD" ]]; then
#    $BRUTE_FORCE_BIN "$MODE" "$LEN_THRESHOLD" < "$INPUT_FILE" > "$BRUTE_OUT"
#else
#    $BRUTE_FORCE_BIN "$MODE" < "$INPUT_FILE" > "$BRUTE_OUT"
#fi
#
#
#tmpfile1=$(mktemp)
#tail -n +5 "$REPEATS_OUT" | sed '$d' > "$tmpfile1"
#header1=$(head -n 1 "$tmpfile1")
#header2=$(head -n 1 "$BRUTE_OUT")
#
#if [[ "$header1" != "$header2" ]]; then
#    echo "Headers do not match."
#    rm "$tmpfile1"
#    exit 1
#fi
#
#rm -f file1.sorted file2.sorted
#
#tail -n +2 "$tmpfile1" | sort -n -k1,1 -k2,2 -k3,3 -k4,4 > file1.sorted
#tail -n +2 "$BRUTE_OUT" | sort -n -k1,1 -k2,2 -k3,3 -k4,4 > file2.sorted
#
#
#
## echo "$header1"
#diff -u file1.sorted file2.sorted
#
#echo "Success! Correct repeats output."
#
#rm "$tmpfile1" "$REPEATS_OUT" "$BRUTE_OUT" "$ROPEBWT3_OUT" "$ROPEBWT3_OUT.optbwtrl"
#rm file1.sorted file2.sorted 
