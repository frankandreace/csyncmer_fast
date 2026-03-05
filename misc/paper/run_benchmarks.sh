#!/usr/bin/env bash
# Collect syncmer detection throughput benchmarks for paper figures.
# Writes tab-separated results to stdout: dataset method throughput_gbps
set -euo pipefail

# ── Paths ──
BENCHMARK=~/tools/csyncmer_fast/misc/tests/benchmark
BENCH_FASTQ=~/tools/csyncmer_fast/misc/fastq/bench_syncmer_fastq
SIMD_MIN=~/tools/simd-minimizers
CHM13=~/data/chm13v2.0.fa
HIFI=~/data/SRR34765324.20G.fastq

# ── Header ──
echo -e "dataset\tmethod\tthroughput_gbps"

# ── digest (actual benchmark) ──
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
bash "$SCRIPT_DIR/benchmark_digest.sh"

# ============================================================================
# (a) CHM13 — whole-genome FASTA, K=31, S=15
# ============================================================================
echo "=== CHM13 (K=31 S=15) ===" >&2

# --- seqhash syng original (Durbin's reference) ---
chm13_seqhash_out=$("$BENCHMARK" --quick "$CHM13" 31 15 seqhash)
echo "$chm13_seqhash_out" >&2

seqhash_mbs=$(echo "$chm13_seqhash_out" | awk '/^SEQHASH_SYNG_ORIG/ {print $3}')
seqhash_gbps=$(awk "BEGIN {printf \"%.4f\", $seqhash_mbs / 1000}")
echo -e "chm13\tseqhash\t$seqhash_gbps"

# --- rescan scalar + twostack SIMD (csyncmer_fast benchmark) ---
chm13_out=$("$BENCHMARK" --quick "$CHM13" 31 15 32)
echo "$chm13_out" >&2

# Parse canonical variants
rescan_mbs=$(echo "$chm13_out" | awk '/^NTH32_CANON_RESCAN_POS/ {print $3}')
rescan_gbps=$(awk "BEGIN {printf \"%.4f\", $rescan_mbs / 1000}")
echo -e "chm13\trescan\t$rescan_gbps"

twostack_mbs=$(echo "$chm13_out" | awk '/^NTH32_CANON_TWOSTACK_POS/ {print $3}')
twostack_gbps=$(awk "BEGIN {printf \"%.4f\", $twostack_mbs / 1000}")
echo -e "chm13\ttwostack\t$twostack_gbps"

# --- multi-8 SIMD (csyncmer_fast benchmark, batches chromosomes 8 at a time) ---
chm13_multi8_out=$("$BENCHMARK" --quick "$CHM13" 31 15 multi8)
echo "$chm13_multi8_out" >&2

multi8_mbs=$(echo "$chm13_multi8_out" | awk '/^NTH32_MULTI8_CANON_POS/ {print $3}')
multi8_gbps=$(awk "BEGIN {printf \"%.4f\", $multi8_mbs / 1000}")
echo -e "chm13\tmulti-8\t$multi8_gbps"

# --- simd-minimizers (FASTA, canonical) ---
echo "--- simd-minimizers CHM13 ---" >&2
simd_out=$("$SIMD_MIN/target/release/examples/chr19_syncmers" "$CHM13")
echo "$simd_out" >&2

# Parse K=31 S=15 block, Canonical pos line: "  Canonical pos:  NNNN syncmers  NNN.N MB/s"
simd_mbs=$(echo "$simd_out" | awk '
    /^K=31 S=15/ { found=1; next }
    found && /Canonical pos/ { print $(NF-1); exit }
')
simd_gbps=$(awk "BEGIN {printf \"%.4f\", $simd_mbs / 1000}")
echo -e "chm13\tsimd-minimizers\t$simd_gbps"

# ============================================================================
# (b) HiFi reads — FASTQ, k=31 w=1022 (K=1052 S=31)
# ============================================================================
echo "=== HiFi reads (K=1052 S=31) ===" >&2

# --- seqhash (Durbin's reference, single-threaded) ---
seqhash_hifi_out=$("$BENCH_FASTQ" -k 31 -w 1022 -seqhash "$HIFI")
echo "$seqhash_hifi_out" >&2
seqhash_hifi_gbps=$(echo "$seqhash_hifi_out" | awk '/^Throughput:/ {print $2}')
echo -e "hifi\tseqhash\t$seqhash_hifi_gbps"

# --- rescan scalar (bench_syncmer_fastq -rescan) ---
rescan_out=$("$BENCH_FASTQ" -k 31 -w 1022 -rescan "$HIFI")
echo "$rescan_out" >&2
rescan_gbps=$(echo "$rescan_out" | awk '/^Throughput:/ {print $2}')
echo -e "hifi\trescan\t$rescan_gbps"

# --- twostack SIMD single (bench_syncmer_fastq -single) ---
twostack_out=$("$BENCH_FASTQ" -k 31 -w 1022 -single "$HIFI")
echo "$twostack_out" >&2
twostack_gbps=$(echo "$twostack_out" | awk '/^Throughput:/ {print $2}')
echo -e "hifi\ttwostack\t$twostack_gbps"

# --- multi-8 split (bench_syncmer_fastq -twopass-nostrand = separate hash + twostack passes) ---
multi8_out=$("$BENCH_FASTQ" -k 31 -w 1022 -twopass-nostrand "$HIFI")
echo "$multi8_out" >&2
multi8_gbps=$(echo "$multi8_out" | awk '/^Throughput:/ {print $2}')
echo -e "hifi\tmulti-8\t$multi8_gbps"

# --- simd-minimizers (FASTQ, canonical) ---
echo "--- simd-minimizers HiFi ---" >&2

# Build bench_syncmer_fastq if needed
FASTQ_EXAMPLE="$SIMD_MIN/examples/bench_syncmer_fastq.rs"
if [ ! -f "$FASTQ_EXAMPLE" ]; then
    echo "Copying bench_syncmer_fastq.rs to simd-minimizers examples..." >&2
    cp ~/tools/csyncmer_fast/misc/fastq/bench_syncmer_fastq.rs "$FASTQ_EXAMPLE"
fi

# Add memmap2 dependency if not present
if ! grep -q 'memmap2' "$SIMD_MIN/Cargo.toml"; then
    echo "Adding memmap2 dependency to simd-minimizers..." >&2
    # Add memmap2 after the [dependencies] section's last entry
    sed -i '/^\[dependencies\]/,/^\[/{/^wide/a memmap2 = "0.9"
    }' "$SIMD_MIN/Cargo.toml"
fi

# Build
echo "Building bench_syncmer_fastq..." >&2
(cd "$SIMD_MIN" && RUSTFLAGS="-C target-cpu=native" cargo build --release --example bench_syncmer_fastq 2>&1) >&2

simd_fastq_out=$("$SIMD_MIN/target/release/examples/bench_syncmer_fastq" -k 31 -w 1022 "$HIFI")
echo "$simd_fastq_out" >&2

# Parse "Throughput: X.XX GB/s"
simd_fastq_gbps=$(echo "$simd_fastq_out" | awk '/^Throughput:/ {print $2}')
echo -e "hifi\tsimd-minimizers\t$simd_fastq_gbps"

echo "=== Done ===" >&2
