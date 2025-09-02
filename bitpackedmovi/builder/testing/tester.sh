#!/bin/bash

if [ $# -ne 3 ]; then 
    echo "please specify exactly 3 command line arguments: builder1, builder2, and generator"
    exit 1
fi

echo $3 $2 $1
