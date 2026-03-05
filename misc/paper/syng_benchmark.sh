#!/usr/bin/env bash
# Benchmark syng end-to-end on HiFi FASTQ: csyncmer vs seqhash (both AVX2=1)
# Usage: ./syng_benchmark.sh [RUNS]
set -euo pipefail

SYNG_DIR=~/tools/syng
HIFI=~/data/SRR34765324.20G.fastq
OUT_DIR=/tmp/syng_hifi_bench
RUNS=${1:-3}

mkdir -p "$OUT_DIR"

echo "========================================"
echo "Syng HiFi Benchmark (multi8 vs csyncmer vs seqhash)"
echo "========================================"
echo "Input: $(basename $HIFI) ($(du -h "$HIFI" | cut -f1))"
echo "Params: -w 1022 -k 31 -T 10"
echo "Runs: $RUNS"
echo ""

# Build all binaries
echo "Building multi8 (AVX2=1 MULTI8=1)..."
make -C "$SYNG_DIR" clean >/dev/null 2>&1
make -C "$SYNG_DIR" AVX2=1 MULTI8=1 -j$(nproc) >/dev/null 2>&1
cp "$SYNG_DIR/syng" "$OUT_DIR/syng_multi8"

echo "Building csyncmer (AVX2=1 CSYNCMER=1)..."
make -C "$SYNG_DIR" clean >/dev/null 2>&1
make -C "$SYNG_DIR" AVX2=1 CSYNCMER=1 -j$(nproc) >/dev/null 2>&1
cp "$SYNG_DIR/syng" "$OUT_DIR/syng_csyncmer"

echo "Building seqhash (AVX2=1)..."
make -C "$SYNG_DIR" clean >/dev/null 2>&1
make -C "$SYNG_DIR" AVX2=1 -j$(nproc) >/dev/null 2>&1
cp "$SYNG_DIR/syng" "$OUT_DIR/syng_seqhash"

echo ""

run_bench() {
    local label=$1
    local bin=$2
    local prefix=$3

    echo "--- $label ($RUNS runs) ---"

    local times=()
    local max_rss=0
    for i in $(seq 1 "$RUNS"); do
        rm -f "$OUT_DIR/${prefix}."*
        /usr/bin/time -f "TIME_OUTPUT %e %M" \
            "$bin" -w 1022 -k 31 -T 10 -o "$OUT_DIR/$prefix" "$HIFI" \
            > "$OUT_DIR/stdout.txt" 2> "$OUT_DIR/time_stderr.txt"
        t=$(grep '^TIME_OUTPUT' "$OUT_DIR/time_stderr.txt" | awk '{print $2}')
        rss=$(grep '^TIME_OUTPUT' "$OUT_DIR/time_stderr.txt" | awk '{print $3}')
        times+=("$t")
        if [ "$rss" -gt "$max_rss" ] 2>/dev/null; then max_rss=$rss; fi
        echo "  Run $i: ${t}s  (RSS: $((rss / 1024))MB)"
    done

    # Sort and pick median
    local sorted=($(printf '%s\n' "${times[@]}" | sort -g))
    local mid=$(( (RUNS - 1) / 2 ))
    echo "  Median: ${sorted[$mid]}s  (min: ${sorted[0]}s, max: ${sorted[$((RUNS-1))]}s)  Peak RSS: $((max_rss / 1024))MB"

    # Syncmer count from last run
    local count
    count=$(grep -oP 'instances of \K\d+(?= syncmers)' "$OUT_DIR/stdout.txt" || echo "?")
    echo "  Syncmers: $count"
    echo ""
}

run_bench "multi8 (ntHash + fastq multi-8)"   "$OUT_DIR/syng_multi8"   "multi8"
run_bench "csyncmer (ntHash + twostack SIMD)" "$OUT_DIR/syng_csyncmer" "csyncmer"
run_bench "seqhash (original hash)"           "$OUT_DIR/syng_seqhash"  "seqhash"

# Cleanup
rm -rf "$OUT_DIR"

echo "Done."
