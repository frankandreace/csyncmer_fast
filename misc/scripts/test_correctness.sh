#!/bin/bash

set -e

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

TEST="../tests/test"
TESTFILE="../tests/test_100kbp.fasta"

# Test with different K and S combinations
for K in 21 31 51 71; do
    for S in 7 11 15 21; do
        if (( S < K )); then
            echo "Test: K=$K, S=$S"
            $TEST "$TESTFILE" $K $S
            echo "---"
        fi
    done
done

echo "All correctness tests passed."
