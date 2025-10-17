#!/bin/bash

set -x
set -e

if [ $# -ne 8 ]; then 
    echo "please specify exactly 8 command line arguments: ropebwt3, BWTbuilder1, BWTbuilder2(kactl), lcpComputer1, lcpComputer2(kactl), generator, loader, and syngenparams"
    exit 1
fi

echo "Making directories for temporary files"
mkdir -p generated
mkdir -p temp
mkdir -p computedOutputs
mkdir -p consoleOutputs


echo "Generating input strings"
$6 $8> generated/stringsAndReverseRaw
head -n -1 generated/stringsAndReverseRaw > generated/stringsAndReverse
awk 'NR%2 == 1' generated/stringsAndReverse > generated/strings    

echo "Generating supporting indexes (ropebwt3 and optbwtrl)"
$1 build -L -d -o temp/a generated/strings > consoleOutputs/ropebwt3
#$2 temp/a temp/a > consoleOutputs/builder

#echo "Generating BWTs"
#$7 BWT temp/a.optbwtrl > computedOutputs/rawBWT1
#$3 < generated/stringsAndReverse > computedOutputs/rawBWT2

#echo "Checking if BWTs match"
#diff computedOutputs/rawBWT1 computedOutputs/rawBWT2

echo "Generating (offset, minLCP) pairs for every run by each method"
#$7 MINLCP temp/a.optbwtrl > computedOutputs/minlcpBuilder
#$4 temp/a computedOutputs/minlcpLCPcomputer > consoleOutputs/lcpComputer
$4 -f fmd -i temp/a -orlcp computedOutputs/minlcpLCPcomputer.rlcp -t temp/tempFile -v quiet
$5 < generated/stringsAndReverse > computedOutputs/minlcpKactl

echo "Checking if minLCP pairs match" 
diff computedOutputs/minlcpKactl computedOutputs/minlcpLCPcomputer.rlcp
#diff computedOutputs/minlcpBuilder computedOutputs/minlcpKactl

echo "Cleaning up files"
rm generated/stringsAndReverseRaw generated/stringsAndReverse generated/strings
rm temp/a consoleOutputs/ropebwt3 #temp/a.optbwtrl temp/a_StructTree.html
#rm computedOutputs/rawBWT1 computedOutputs/rawBWT2
#rm computedOutputs/minlcpLCPcomputer.optbwtrl computedOutputs/minlcpLCPcomputer_StructTree.html 
rm computedOutputs/minlcpLCPcomputer.rlcp computedOutputs/minlcpKactl #computedOutputs/minlcpBuilder 
rm construction.html #consoleOutputs/builder consoleOutputs/lcpComputer 
