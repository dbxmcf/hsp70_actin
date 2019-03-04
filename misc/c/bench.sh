#!/bin/bash
arr=( 1 2 5 10 20 )
for i in "${arr[@]}";
do
    export OMP_NUM_THREADS="$i"
    SECONDS=0
    time ./a.out
    echo "threads=$OMP_NUM_THREADS, $SECONDS sec"
done
