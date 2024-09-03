#!/usr/bin/env bash
rm Results/*dM.txt
ncores=10

for i in $(seq 1 $ncores)
do
    let j="$i"-1
    dir="Part_$j"
    cd $dir
    cat ExpanddM.txt >> ../Results/ExpanddM.txt
    cat TestdM.txt >> ../Results/BorndM.txt
    cat ResumdM.txt >> ../Results/ResumdM.txt
    
    cd ..
    ##rm -r $dir
done
