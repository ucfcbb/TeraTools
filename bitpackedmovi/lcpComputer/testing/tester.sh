#!/bin/bash

set -x
set -e

if [ $# -ne 5 ]; then 
    echo "please specify exactly 5 command line arguments: ropebwt3, lcpComputer1, lcpComputer2(kactl), generator, and syngenparams"
    exit 1
fi

ropebwt3=$1
lcpComputer1=$2
lcpComputer2_kactl=$3
generator=$4
syngenparams=$5

echo "Making directories for temporary files"
mkdir -p generated
mkdir -p temp
mkdir -p computedOutputs
mkdir -p consoleOutputs


echo "Generating input strings"
$generator $syngenparams> generated/stringsAndReverseRaw
head -n -1 generated/stringsAndReverseRaw > generated/stringsAndReverse
awk 'NR%2 == 1' generated/stringsAndReverse > generated/strings    

echo "Generating supporting index (ropebwt3)"
$ropebwt3 build -L -d -o temp/a generated/strings > consoleOutputs/ropebwt3

echo "Generating (offset, minLCP) pairs for every run by each method"
$lcpComputer1 -f fmd -i temp/a -orlcp computedOutputs/minlcpLCPcomputer -t temp/tempFile -v quiet
$lcpComputer2_kactl < generated/stringsAndReverse > computedOutputs/minlcpKactl

echo "Checking if minLCP pairs match" 
diff computedOutputs/minlcpKactl computedOutputs/minlcpLCPcomputer.rlcp

echo "Cleaning up files"
rm generated/stringsAndReverseRaw generated/stringsAndReverse generated/strings
rm temp/a consoleOutputs/ropebwt3
rm computedOutputs/minlcpLCPcomputer.rlcp computedOutputs/minlcpKactl
rm construction.html
