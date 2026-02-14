#!/bin/bash

set -e

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

# Default values
DEFAULT_FILE="../tests/test_100kbp.fasta"
DEFAULT_KMER_SIZE=31
DEFAULT_SMER_SIZE=11
CLUSTER_MODE=0

# Parse command-line options
while getopts "f:k:s:c" opt; do
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
    c)
      CLUSTER_MODE=1
      ;;
    \?)
      echo "Invalid option: -$OPTARG" >&2
      echo "Usage: $0 [-f fasta_file] [-k kmer_size] [-s smer_size] [-c]"
      echo "  -c  Cluster mode: skip CPU frequency/SMT manipulation (no sudo)"
      exit 1
      ;;
  esac
done

FILE="${FILE:-$DEFAULT_FILE}"
KMER_SIZE="${KMER_SIZE:-$DEFAULT_KMER_SIZE}"
SMER_SIZE="${SMER_SIZE:-$DEFAULT_SMER_SIZE}"
OUTDIR="../tests/results"
OUTFILE="$OUTDIR/benchmark.tsv"

BENCHMARK="../tests/benchmark"

echo "TESTING SPEED"
mkdir -p "$OUTDIR"

if [ "$CLUSTER_MODE" -eq 0 ]; then
    echo "Setting cpu to 2.6 GHz."
    echo "[Executing] sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz"
    sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz

    SMT_STATUS=$(cat /sys/devices/system/cpu/smt/control)

    if [ "$SMT_STATUS" = "on" ]; then
        echo "Disabling HYPERTHREADING for more accurate time estimation"
        echo "[Executing] sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'"
        sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'
    else
        echo "HYPERTHREADING is already disabled (status: $SMT_STATUS)."
    fi
fi

echo "RUNNING SPEED TEST"
echo "[Executing] $BENCHMARK $FILE $KMER_SIZE $SMER_SIZE $OUTFILE"
$BENCHMARK $FILE $KMER_SIZE $SMER_SIZE $OUTFILE

if [ "$CLUSTER_MODE" -eq 0 ] && [ "$SMT_STATUS" = "on" ]; then
    echo "Re-enabling HYPERTHREADING."
    echo "[Executing] sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'"
    sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'
fi
