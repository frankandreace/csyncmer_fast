#!/bin/bash

set -e 


# Default values
DEFAULT_FILE="data/chr19_bit.fa"
DEFAULT_KMER_SIZE=31
DEFAULT_SMER_SIZE=11
DEFAULT_MODE=1

# Parse command-line options
while getopts "f:k:s:r:" opt; do
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
    r)
      MODE="$OPTARG"
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
MODE="${MODE:-$DEFAULT_MODE}"

echo "TESTING SPEED"
mkdir -p ./benchmark

echo "RUNNING SPEED TEST"
# echo "[Executing] ./bin/test $FILE $KMER_SIZE $SMER_SIZE $MODE"
./bin/test $FILE $KMER_SIZE $SMER_SIZE $MODE
