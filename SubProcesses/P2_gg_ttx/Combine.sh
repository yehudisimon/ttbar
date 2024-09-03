#!/usr/bin/env bash
rm Results/*dM.txt
#rm Results/*.txt
ncores=10

for i in $(seq 1 $ncores)
do
    let j="$i"-1
    dir="Part_$j"
    cd $dir
    cat ExpanddM_gg.txt >> ../Results/ExpanddM.txt
    cat ResumdM_gg.txt >> ../Results/ResumdM.txt
    cat BorndM_gg.txt >> ../Results/BorndM.txt
    cat DiffdM_gg.txt >> ../Results/DiffdM.txt
    # cat BorndM_gg.txt >> ../Results/Born_alpha.txt
    # cat Diff_gg.txt >> ../Results/Diff_alpha.txt
    # cat Expand_gg.txt >> ../Results/Expand_alpha.txt
    
    
    cd ..
    ##rm -r $dir
done
