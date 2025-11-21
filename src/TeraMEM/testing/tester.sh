#!/usr/bin/env bash


#set -x
set -e

if [[ $# -lt 9 ]]; then
    echo "please specify exactly 9 command line arguments: ropebwt3, TeraLCP, msIndexBuilder, MEM, bruteForceRepeats, generator, syngenparams, mode:[SM,LM], lengthThreshold"
    exit 1
fi

ropebwt3=$1
TeraLCP=$2
msIndexBuilder=$3
MEM=$4
MEM_BRUTE=$5
generator=$6
syngenparams=$7
MODE=$8
LEN=$9

# Validate mode parameter
if [[ "$MODE" != "SM" && "$MODE" != "LM" ]]; then
    echo "Error: Mode must be either 'SM' or 'LM'"
    exit 1
fi

echo "Making directories for temporary files"
mkdir -p generated
mkdir -p temp
mkdir -p computedOutputs
mkdir -p consoleOutputs

echo "Generating input strings"
$generator $syngenparams > generated/stringsAndReverseRaw
head -n -1 generated/stringsAndReverseRaw > generated/stringsAndReverse
awk 'NR%2 == 1' generated/stringsAndReverse > generated/strings    

echo "Generating supporting index (ropebwt3)"
$ropebwt3 build -L -d -o temp/a.fmd generated/strings

MEM_OUT="computedOutputs/MEM.out"
BRUTE_OUT="computedOutputs/bruteMEM.out"

echo "Running TeraLCP on temp/a.fmd"
$TeraLCP -f fmd -i temp/a.fmd -t temp/temp -oindex temp/a -v quiet

echo "Running msIndexBuilder on temp/a.lcp_index"
$msIndexBuilder -i temp/a.lcp_index -o temp/a -v quiet

echo "Running MEM on temp/a.ms_index"
$MEM -MEM $MODE -i "temp/a.ms_index" -L $LEN -o $MEM_OUT -v quiet

echo "Running brute_force_repeats on generated/stringsAndReverse"
$MEM_BRUTE $MODE $LEN < generated/stringsAndReverse > $BRUTE_OUT


header1=$(head -n 1 "$MEM_OUT")
header2=$(head -n 1 "$BRUTE_OUT")

if [[ "$header1" != "$header2" ]]; then
    echo "Headers do not match."
    rm "$tmpfile1"
    exit 1
fi

sort $MEM_OUT   -n -k1,1 -k2,2 -k3,3 -k4,4 > $MEM_OUT.sorted
sort $BRUTE_OUT -n -k1,1 -k2,2 -k3,3 -k4,4 > $BRUTE_OUT.sorted

diff -u $MEM_OUT.sorted $BRUTE_OUT.sorted

#echo "Success! Correct repeats output."

echo "Cleaning up files"
rm generated/stringsAndReverseRaw generated/stringsAndReverse generated/strings
rm temp/a.fmd temp/a.lcp_index temp/temp construction.html temp/a.ms_index
rm $MEM_OUT $BRUTE_OUT
rm $MEM_OUT.sorted $BRUTE_OUT.sorted 
