#!/bin/bash

thread_nums=(1 2 4 8 16)

for i in ${thread_nums[@]}; do
    export OMP_NUM_THREADS=$i
    ./hot
done
