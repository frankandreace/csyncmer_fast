#!/bin/bash

set -e

# GETTING THE DIRECTORY OF THE SCRIPT
SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

# TESTING RESULTS ARE CORRECT
# ./test_correctness.sh


# TESTING COMPUTATION SPEED
OUTDIR="../tests/results"
rm -rf "$OUTDIR/benchmark.tsv" "$OUTDIR/benchmark_plot"*

for i in {1..10}
do
    ./test_speed.sh
done

./plot_result.py "$OUTDIR/benchmark.tsv" "$OUTDIR/benchmark_plot"
