#!/bin/bash

set -e

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

# Default values
DEFAULT_FILE="../../data/chr1.fa"
DEFAULT_KMER_SIZE=31
DEFAULT_SMER_SIZE=11

# Parse command-line options
while getopts "f:k:s:" opt; do
  case $opt in
    f)
      FILE="$OPTARG"
      ;;
    k)
      KMER_SIZE="$OPTARG"
      ;;
    s)
      SMER_SIZE="$OPTARG"
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      exit 1
      ;;
  esac
done

FILE="${FILE:-$DEFAULT_FILE}"
KMER_SIZE="${KMER_SIZE:-$DEFAULT_KMER_SIZE}"
SMER_SIZE="${SMER_SIZE:-$DEFAULT_SMER_SIZE}"
OUTFILE="../../benchmark/results/benchmark.tsv"

BENCHMARK="../bench/benchmark"

echo "TESTING SPEED"
mkdir -p "../../benchmark/results"

echo "RUNNING SPEED TEST"
$BENCHMARK bench $FILE $KMER_SIZE $SMER_SIZE $OUTFILE
