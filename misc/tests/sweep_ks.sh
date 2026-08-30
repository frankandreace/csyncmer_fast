#!/usr/bin/env bash
# Sweep syncmer detection throughput over a range of commonly used (k, s) values.
# Writes TSV to stdout: k s w method syncmers throughput_gbps
set -euo pipefail

BENCHMARK=~/tools/csyncmer_fast/misc/tests/benchmark
SIMD_MIN=~/tools/simd-minimizers/target/release/examples/sweep_syncmers
FASTA=${1:-~/data/chm13v2.0.fa}

# Commonly used (k, s) pairs; w = k - s + 1
CONFIGS="15,7 15,11 21,11 21,15 31,15 31,19 31,23 41,21 51,31"

echo -e "k\ts\tw\tmethod\tsyncmers\tthroughput_gbps"

for cfg in $CONFIGS; do
    k=${cfg%,*}; s=${cfg#*,}; w=$((k - s + 1))
    echo "=== k=$k s=$s w=$w ===" >&2

    out32=$("$BENCHMARK" --quick "$FASTA" "$k" "$s" 32)
    echo "$out32" >&2
    echo "$out32" | awk -v k=$k -v s=$s -v w=$w '
        /^NTH32_CANON_RESCAN_POS/   {printf "%s\t%s\t%s\trescan\t%s\t%.4f\n",   k,s,w,$2,$3/1000}
        /^NTH32_CANON_TWOSTACK_POS/ {printf "%s\t%s\t%s\ttwostack\t%s\t%.4f\n", k,s,w,$2,$3/1000}'

    out8=$("$BENCHMARK" --quick "$FASTA" "$k" "$s" multi8)
    echo "$out8" >&2
    echo "$out8" | awk -v k=$k -v s=$s -v w=$w '
        /^NTH32_MULTI8_CANON_POS/ {printf "%s\t%s\t%s\tmulti-8\t%s\t%.4f\n", k,s,w,$2,$3/1000}'

    outsm=$("$SIMD_MIN" "$FASTA" "$k" "$s")
    echo "$outsm" >&2
    echo "$outsm" | awk -v k=$k -v s=$s -v w=$w '
        /Canonical pos/ {printf "%s\t%s\t%s\tsimd-minimizers\t%s\t%.4f\n", k,s,w,$3,$(NF-1)/1000}'
done

echo "=== sweep done ===" >&2
