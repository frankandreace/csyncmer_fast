#!/bin/bash

set -e 

# GETTING THE DIRECTORY OF THE SCRIPT
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

#From now on all the paths are relative to here

# TESTING RESULTS ARE CORRECT
# ./scripts/test_correctness.sh


# TESTING COMPUTATION SPEED
rm -rf benchmark/results/benchmark.tsv benchmark/results/benchmark_plot*

for i in {1..10} 
do
    ./scripts/test_speed.sh
done
# ./scripts/test_speed.sh
./scripts/plot_result.py benchmark/results/benchmark.tsv benchmark/results/benchmark_plot