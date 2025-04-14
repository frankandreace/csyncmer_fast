#!/bin/bash

set -e 

# TESTING RESULTS ARE CORRECT
# ./test_correctness.sh


# TESTING COMPUTATION SPEED
rm -rf benchmark/benchmark.tsv benchmark/benchmark_plot*

for i in {1..20} 
do
    ./scripts/test_speed.sh
done

./scripts/plot_result.py benchmark/benchmark.tsv benchmark/benchmark_plot