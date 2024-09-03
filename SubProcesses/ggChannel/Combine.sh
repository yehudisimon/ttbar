#!/usr/bin/env bash
rm Results/*.txt
ncores=10

for i in $(seq 1 $ncores)
do
    let j="$i"-1
    dir="Part_$j"
    cd $dir
    cat BorndM_gg.txt >> ../Results/BorndM.txt
    cat ResumdM_gg.txt >> ../Results/ResumdM.txt
    
    cd ..
done
