#!/bin/bash

set -e 

# TESTING RESULTS ARE CORRECT
./scripts/test_correctness.sh


# TESTING COMPUTATION SPEED

mkdir -p ./benchmark
rm -rf benchmark/benchmark.tsv benchmark/benchmark_plot*

for i in {1..10} 
do
    ./scripts/test_speed.sh
done

./scripts/plot_result.py benchmark/benchmark.tsv benchmark/benchmark_plot