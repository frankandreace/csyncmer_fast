#!/usr/bin/env bash
# Sweep syncmer detection throughput over (k, s) values on FASTQ long reads.
# Writes TSV to stdout: k s w method syncmers throughput_gbps
set -euo pipefail

BENCH=~/tools/csyncmer_fast/misc/fastq/bench_syncmer_fastq
SIMD_MIN=~/tools/simd-minimizers/target/release/examples/bench_syncmer_fastq
FASTQ=${1:-~/data/SRR34765324.20G.fastq}

# Mirrors the CHM13 grid (K,s), plus the large-window configs used for long reads.
# w = K - s + 1; bench takes -k <s> -w <w>.
CONFIGS="15,7 15,11 21,11 21,15 31,15 31,19 31,23 41,21 51,31 285,31 541,31 1052,31"

echo -e "k\ts\tw\tmethod\tsyncmers\tthroughput_gbps"

parse() {  # stdin -> "syncmers throughput"
    awk '/^Syncmers:/ {c=$2} /^Throughput:/ {t=$2} END {print c, t}'
}

emit() {  # K s w method <<< output
    local K=$1 s=$2 w=$3 m=$4
    read -r cnt tp < <(parse)
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" "$K" "$s" "$w" "$m" "$cnt" "$tp"
}

for cfg in $CONFIGS; do
    K=${cfg%,*}; s=${cfg#*,}; w=$((K - s + 1))
    echo "=== K=$K s=$s w=$w ===" >&2

    for mode in rescan:rescan single:twostack twopass-nostrand:multi-8; do
        flag=${mode%:*}; name=${mode#*:}
        out=$("$BENCH" -k "$s" -w "$w" "-$flag" "$FASTQ")
        echo "$out" >&2
        echo "$out" | emit "$K" "$s" "$w" "$name"
    done

    out=$("$SIMD_MIN" -k "$s" -w "$w" "$FASTQ")
    echo "$out" >&2
    echo "$out" | emit "$K" "$s" "$w" "simd-minimizers"
done

echo "=== sweep done ===" >&2
