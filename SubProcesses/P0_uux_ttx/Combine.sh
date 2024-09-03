#!/usr/bin/env bash
rm Results/*dM.txt
#rm Results/*.txt
ncores=10

for i in $(seq 1 $ncores)
do
    let j="$i"-1
    dir="Part_$j"
    cd $dir
    cat ExpanddM_qqb.txt >> ../Results/ExpanddM.txt
    cat ResumdM_qqb.txt >> ../Results/ResumdM.txt
    cat BorndM_qqb.txt >> ../Results/BorndM.txt
    cat DiffdM_qqb.txt >> ../Results/DiffdM.txt
    # cat BorndM_qqb.txt >> ../Results/Born_alpha.txt
    # cat Diff_qqb.txt >> ../Results/Diff_alpha.txt
    
    
    cd ..
    ##rm -r $dir
done
