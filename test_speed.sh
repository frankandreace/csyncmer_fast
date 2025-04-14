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

echo "Setting cpu to 2.6 GHz."
echo "[Executing] sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz"
sudo cpupower frequency-set --governor powersave -d 2.6GHz -u 2.6GHz

echo "Disabling HYPERTREADING for more accurate time estimation"
echo "[Executing] sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'"
sudo sh -c 'echo off > /sys/devices/system/cpu/smt/control'

echo "RUNNING SPEED TEST"
echo "[Executing] $FILE $KMER_SIZE $SMER_SIZE $MODE"
./test $FILE $KMER_SIZE $SMER_SIZE $MODE

echo "Re-enabling HYPERTREADING."
echo "[Executing] sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'"
sudo sh -c 'echo on > /sys/devices/system/cpu/smt/control'