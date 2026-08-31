#!/usr/bin/env bash
# Benchmark syng end-to-end on HPRC HiFi reads: multi8 vs csyncmer vs seqhash.
# Produces the data behind panel (c) of the paper figure, in syng.tsv's schema.
#
# Recorded run (the numbers currently in syng.tsv):
#   10 DeepConsensus HiFi read files from HPRC Release 2 -- one PacBio Sequel II
#   SMRT cell from each of 10 samples, 636 GB / 317.7 Gbp total
#   syncmer length 1023, s-mer 31, 8 threads, 1 run
#   AMD EPYC 7552 node
#     multi8    581.83s   315006298 syncmers
#     csyncmer  597.05s   315006298 syncmers
#     seqhash  1177.57s   328409682 syncmers
#
# Usage: ./syng_benchmark.sh [RUNS] [THREADS]
set -euo pipefail

SYNG_DIR=~/tools/syng
DATA_DIR=${DATA_DIR:-~/data/hprc}
OUT_DIR=${OUT_DIR:-/tmp/syng_hifi_bench}
RUNS=${1:-1}
THREADS=${2:-8}

# ── Parameters ──
# syng and csyncmer_fast's bench_syncmer_fastq need -w values ONE APART for the
# same syncmer length. seqhash's internal w counts s-mer hashes, so the syncmer
# length is always sh->w + sh->k - 1; syng calls seqhashCreate(k, w+1) with an
# extra +1 (syng.c:542, and syngBWTcreate uses k+w at syng.c:509), giving
# syncmer length = w + k, while bench_syncmer_fastq passes the internal window
# straight through, giving K = w + k - 1. So syncmer length 1023 at s=31 is
# syng -w 992 and bench_syncmer_fastq -w 993. Keep them in step.
SMER=31
SYNCMER_LEN=1023
W=$((SYNCMER_LEN - SMER))

# ── Input: one DeepConsensus HiFi file per HPRC Release 2 sample ──
# Index: human-pangenomics/hprc_intermediate_assembly,
#        data_tables/sequencing_data/data_hifi_release2_v1.0.index.csv
# (each of these samples has 3-6 DeepConsensus files; this is one cell each)
SAMPLES=(
    HG00099.m54329U_220825_174247.dc.q20.fastq
    HG00280.m54329U_220901_221341.dc.q20.fastq
    HG00558.m54329U_220107_233847.dc.q20.fastq
    HG00639.m54329U_211222_104516.dc.q20.fastq
    HG01074.m54329U_211110_112322.dc.q20.fastq
    HG01123.m54329U_200205_002609.dc.q20.fastq
    HG02257.m64076_200125_231256.dc.q20.fastq
    HG02486.m64076_200211_192227.dc.q20.fastq
    HG02922.m54329U_220816_182601.dc.q20.fastq
    HG03209.m64076_220526_115049.dc.q20.fastq
)

INPUTS=()
for s in "${SAMPLES[@]}"; do
    f="$DATA_DIR/$s"
    [ -f "$f" ] || { echo "missing input: $f" >&2; exit 1; }
    INPUTS+=("$f")
done

# Pin to CPUs 0..THREADS-1 for reproducibility
CPUS="0-$((THREADS - 1))"
[ "$THREADS" -eq 1 ] && CPUS="0"
TASKSET="taskset -c $CPUS"

mkdir -p "$OUT_DIR"

total_bytes=$(du -cb "${INPUTS[@]}" | tail -1 | cut -f1)

echo "========================================"
echo "Syng HPRC HiFi Benchmark (multi8 vs csyncmer vs seqhash)"
echo "========================================"
echo "Input:  ${#INPUTS[@]} files, $((total_bytes / 1000000000)) GB"
echo "Params: -w $W -k $SMER (syncmer length $SYNCMER_LEN) -T $THREADS (pinned to CPUs $CPUS)"
echo "Runs:   $RUNS"
echo ""

# Build all binaries
echo "Building multi8 (AVX2=1 MULTI8=1)..."
make -C "$SYNG_DIR" clean >/dev/null 2>&1
make -C "$SYNG_DIR" AVX2=1 MULTI8=1 -j"$(nproc)" >/dev/null 2>&1
cp "$SYNG_DIR/syng" "$OUT_DIR/syng_multi8"

echo "Building csyncmer (AVX2=1 CSYNCMER=1)..."
make -C "$SYNG_DIR" clean >/dev/null 2>&1
make -C "$SYNG_DIR" AVX2=1 CSYNCMER=1 -j"$(nproc)" >/dev/null 2>&1
cp "$SYNG_DIR/syng" "$OUT_DIR/syng_csyncmer"

echo "Building seqhash (AVX2=1)..."
make -C "$SYNG_DIR" clean >/dev/null 2>&1
make -C "$SYNG_DIR" AVX2=1 -j"$(nproc)" >/dev/null 2>&1
cp "$SYNG_DIR/syng" "$OUT_DIR/syng_seqhash"

echo ""

# config → "median min max syncmers", filled in by run_bench
declare -A RESULT

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
            $TASKSET "$bin" -w "$W" -k "$SMER" -T "$THREADS" \
                -o "$OUT_DIR/$prefix" "${INPUTS[@]}" \
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

    RESULT[$prefix]="${sorted[$mid]} ${sorted[0]} ${sorted[$((RUNS-1))]} $count"
}

run_bench "multi8 (ntHash + fastq multi-8)"   "$OUT_DIR/syng_multi8"   "multi8"
run_bench "csyncmer (ntHash + twostack SIMD)" "$OUT_DIR/syng_csyncmer" "csyncmer"
run_bench "seqhash (original hash)"           "$OUT_DIR/syng_seqhash"  "seqhash"

# ── Emit syng.tsv (the input to plot_throughput.py -s) ──
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
TSV="$SCRIPT_DIR/syng.tsv"
read -r base_median _ _ _ <<< "${RESULT[seqhash]}"
{
    printf "config\tmedian_s\tmin_s\tmax_s\tvs_seqhash\tsyncmers\tthreads\tpinned\n"
    for cfg in multi8 csyncmer seqhash; do
        read -r med mn mx cnt <<< "${RESULT[$cfg]}"
        if [ "$cfg" = seqhash ]; then
            vs=baseline
        else
            vs=$(awk "BEGIN {printf \"%+.0f%%\", 100 * ($med - $base_median) / $base_median}")
        fi
        printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$cfg" "$med" "$mn" "$mx" "$vs" "$cnt" "$THREADS" "$CPUS"
    done
} > "$TSV.tmp" && mv "$TSV.tmp" "$TSV"
echo "Wrote $TSV"

rm -rf "$OUT_DIR"
echo "Done."
