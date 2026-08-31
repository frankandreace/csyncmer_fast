#!/usr/bin/env bash
# Benchmark digest library syncmer throughput on CHM13 and HiFi datasets.
# Outputs TSV lines: dataset method throughput_gbps
set -euo pipefail

DIGEST_DIR=~/tools/digest
CHM13=~/data/chm13v2.0.fa
HIFI=~/data/SRR34765324.20G.fastq
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BENCH="$SCRIPT_DIR/benchmark_digest"
RUNS=3
SECTION="${SECTION:-all}"

# ── Build ──
if [ ! -f "$BENCH" ] || [ "$SCRIPT_DIR/benchmark_digest.cpp" -nt "$BENCH" ]; then
    echo "Building benchmark_digest ..." >&2
    g++ -std=c++17 -O3 -o "$BENCH" "$SCRIPT_DIR/benchmark_digest.cpp" \
        -I "$DIGEST_DIR/build/include" -L "$DIGEST_DIR/build/lib" -lnthash
fi

# ── CHM13: K=31, S=15 → k_small=15, large_window=17 ──
if [ "$SECTION" = all ] || [ "$SECTION" = chm13 ]; then
    echo "=== digest CHM13 (K=31 S=15) ===" >&2
    chm13_gbps=$("$BENCH" "$CHM13" 15 17 "$RUNS")
    echo -e "chm13\tdigest\t$chm13_gbps"
fi

# ── HiFi: K=1023, S=31 → k_small=31, large_window=993 ──
if [ "$SECTION" = all ] || [ "$SECTION" = hifi ]; then
    echo "=== digest HiFi (K=1023 S=31) ===" >&2
    hifi_gbps=$("$BENCH" "$HIFI" 31 993 "$RUNS")
    echo -e "hifi\tdigest\t$hifi_gbps"
fi

echo "=== digest done ===" >&2
