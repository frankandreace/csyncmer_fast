#!/bin/bash

set -e

# GETTING THE DIRECTORY OF THE SCRIPT
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

# TESTING RESULTS ARE CORRECT
# ./test_correctness.sh


# TESTING COMPUTATION SPEED
rm -rf ../../benchmark/results/benchmark.tsv ../../benchmark/results/benchmark_plot*

for i in {1..10}
do
    ./test_speed.sh
done

./plot_result.py ../../benchmark/results/benchmark.tsv ../../benchmark/results/benchmark_plot
