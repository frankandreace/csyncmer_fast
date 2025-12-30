#!/bin/bash

set -e

SCRIPT_DIR=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )

cd "$SCRIPT_DIR"

# Default values
DEFAULT_FILE="../../data/chr19_bit.fa"
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

echo "RUNNING SPEED TEST"
echo "[Executing] $BENCHMARK bench $FILE $KMER_SIZE $SMER_SIZE $OUTFILE"
$BENCHMARK bench $FILE $KMER_SIZE $SMER_SIZE $OUTFILE

if [ "$SMT_STATUS" = "on" ]; then
    echo "Re-enabling HYPERTREADING."
    echo "[Executing] sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'"
    sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'
fi
