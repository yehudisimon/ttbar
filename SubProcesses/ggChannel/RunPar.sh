ncores=10

for i in $(seq 1 $ncores)
do    
    let j="$i"-1
    dir="Part_$j"
    if [ -d "${dir}" ];
    then
    	rm -rf $dir
    fi
    mkdir $dir
    cd $dir
    ln -s ../RUN ./
    nohup ./RUN 21200 $i $ncores &# > output.log &
    echo "PID of process n°" $i " = " $! 
    cd ../
done
