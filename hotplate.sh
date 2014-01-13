#!/bin/bash

echo "Stuff"
export OMP_NUM_THREADS=1
./hot

export OMP_NUM_THREADS=2
./hot

export OMP_NUM_THREADS=4
./hot

export OMP_NUM_THREADS=8
./hot

export OMP_NUM_THREADS=16
./hot
