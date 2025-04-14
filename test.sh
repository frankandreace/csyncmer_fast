#!/bin/bash

set -e 

# TESTING RESULTS ARE CORRECT
# ./test_correctness.sh


# TESTING COMPUTATION SPEED

for i in {1..20} 
do
    ./test_speed.sh
done

./plot_result.py benchmark.tsv benchmark_plot