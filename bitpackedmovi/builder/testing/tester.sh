#!/bin/bash

if [ $# -ne 4 ]; then 
    echo "please specify exactly 3 command line arguments: ropebwt3, builder1, builder2, and generator"
    exit 1
fi

mkdir -p generated
mkdir -p temp
mkdir -p bwt

$4 > generated/stringsAndReverse
awk 'NR%2 == 1' generated/stringsAndReverse > generated/strings    

$3 generated/stringsAndReverse > bwt/a_kactl

$1 build -L -d -o temp/a generated/strings
$2 temp/a > bwt/a_ours

echo $4 $3 $2 $1
