#!/bin/bash

thread_nums=(1 2 4 8 16 32)

for i in ${thread_nums[@]}; do
    ./hot $i 10
done
