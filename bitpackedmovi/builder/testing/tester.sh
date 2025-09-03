#!/bin/bash

set -x
set -e

if [ $# -ne 4 ]; then 
    echo "please specify exactly 3 command line arguments: ropebwt3, builder1, builder2, and generator"
    exit 1
fi

mkdir -p generated
mkdir -p temp
mkdir -p bwt

$4 > generated/stringsAndReverseRaw
tail -n +3 generated/stringsAndReverseRaw > generated/stringsAndReverse
#tail -n +3 generated/stringsAndReverseRaw > generated/stringsAndReverseMetadata
#head -n -1 generated/stringsAndReverseMetadata > generated/stringsAndReverse

awk 'NR%2 == 1' generated/stringsAndReverse > generated/strings    

head -n 2 generated/stringsAndReverseRaw > generated/stringsAndReverseMetadata
#tail -n 1 generated/stringsAndReverseRaw >> generated/stringsAndReverseMetadata

cat <(wc -l < generated/stringsAndReverse) generated/stringsAndReverse | $3

$1 build -L -d -o temp/a generated/strings
$2 temp/a temp/a > bwt/a_ours

echo $4 $3 $2 $1
