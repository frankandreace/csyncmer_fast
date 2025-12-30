#!/bin/bash

set -e

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

BENCHMARK="../bench/benchmark"

generate_sequence() {
    local len=$1
    local seq=""
    local bases=(A C G T)
    for (( i=0; i<$len; i++ )); do
        seq+=${bases[$RANDOM % 4]}
    done
    echo "$seq"
}

a=5
b=2

for i in {1..10}
do
    a=$(( RANDOM % 989 + 11 ))  # Random integer between 11 and 999 inclusive

    # Calculate 'b' constraints: 6 < b < a and b < 100
    b_min=2
    b_max=$(( a - 1 ))
    if (( b_max > 31 )); then
        b_max=31
    fi

    # Ensure b_min does not exceed b_max
    if (( b_max < b_min )); then
        b=$b_min
    else
        b=$(( RANDOM % (b_max - b_min + 1) + b_min ))
    fi
    sequence=$(generate_sequence 10000)
    $BENCHMARK check "$sequence" $a $b
done

echo "All tests for correctness are OK."
