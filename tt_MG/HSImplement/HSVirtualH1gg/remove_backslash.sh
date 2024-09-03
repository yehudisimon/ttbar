#!/usr/bin/env bash

export LANG=C
for file in $(ls hardfuncs*.f)
do
#    echo $file
    sed 's|\\$| |g' $file> temp.f
    mv temp.f $file
done