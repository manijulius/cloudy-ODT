#! /bin/bash

in_dir=$1
nml=$2
exp_name=$3

n_realization=5
n_batch=1
let n_realz_in_batch=$n_realization/$n_batch


for i in $(seq 1 $n_batch)
do
    
    echo "********** Batch :$i "
    for j in $(seq 1 $n_realz_in_batch)
    do
        
        let n=($i-1)*$n_realz_in_batch+$j
        realz_idx=$(printf "%03d" $n)
        echo "Realization : $realz_idx"
        nohup ./LabExp_mphy.exe $in_dir $nml $exp_name $realz_idx &> ${exp_name}_realz_${realz_idx}.log &
        echo "Process id: $!"       
    done;
    
#    sleep 25000 
    echo "waiting for batch $i processes to complete"
    wait;

done;
