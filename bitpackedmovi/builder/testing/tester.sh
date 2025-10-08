#!/bin/bash

set -x
set -e

if [ $# -ne 6 ]; then 
    echo "please specify exactly 6 command line arguments: ropebwt3, builder1, builder2, generator, loader, and parameters for generator"
    exit 1
fi

mkdir -p generated
mkdir -p temp
mkdir -p bwt

$4 $6 > generated/stringsAndReverseRaw
head -n -1 generated/stringsAndReverseRaw > generated/stringsAndReverse
#tail -n +3 generated/stringsAndReverseRaw > generated/stringsAndReverseMetadata
#head -n -1 generated/stringsAndReverseMetadata > generated/stringsAndReverse

awk 'NR%2 == 1' generated/stringsAndReverse > generated/strings    

tail -n 2 generated/stringsAndReverseRaw > generated/stringsAndReverseMetadata
#tail -n 1 generated/stringsAndReverseRaw >> generated/stringsAndReverseMetadata

#cat <(wc -l < generated/stringsAndReverse) generated/stringsAndReverse | $3
$3 < generated/stringsAndReverse > temp/rawBWT2

$1 build -L -d -o temp/a generated/strings
$2 temp/a temp/a > bwt/a_ours

$5 BWT temp/a.optbwtrl > temp/rawBWT1

diff temp/rawBWT1 temp/rawBWT2

# if [ $? -ne 0 ]; then
#     echo "mismatch between linear SA vs. run length compressed SA"
#     exit 1
# fi

#echo $5 $4 $3 $2 $1

rm generated/stringsAndReverseRaw generated/stringsAndReverse generated/strings generated/stringsAndReverseMetadata
rm temp/{a,a.optbwtrl,a_StructTree.html} bwt/a_ours
rm temp/rawBWT{1,2}
