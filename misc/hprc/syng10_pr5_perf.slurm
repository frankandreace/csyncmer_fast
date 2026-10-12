#!/usr/bin/env bash
#SBATCH -J syng10_pr5
#SBATCH -p seqbio
#SBATCH -x maestro-2499
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH -t 12:00:00
#SBATCH -o /pasteur/appa/scratch/rchikhi/hprc_v2/reads/logs/syng10pr5_%j.out
# Fig 1(c) from syng itself: PR5 (rchikhi/opt-avx2-csyncmer @ 698d409) end-to-end on the
# 10 HPRC HiFi files, -w 992 -k 31 -T 8 pinned to CPUs 0-7. Per config: one plain timed
# run (the bar total) and one perf run (extraction share of samples, for the split).
# PR5 has no separate multi8 build: its CSYNCMER path is the 8-lane multi-read one.
set -euo pipefail
SYNG=/pasteur/helix/users/rchikhi/tools/syng
COMMIT=698d409
SRC=/pasteur/appa/scratch/rchikhi/hprc_v2/reads
OUT=$SRC/logs/syng10pr5_$SLURM_JOB_ID
W=/dev/shm/syng10pr5_$SLURM_JOB_ID
trap 'rm -rf "$W"' EXIT
mkdir -p "$W" "$OUT"

echo "== node $(hostname) $(date)"; lscpu | grep -E "Model name|Socket"
git -C "$SYNG" log --oneline -1 "$COMMIT"

MANIFEST=/pasteur/helix/users/rchikhi/data/hprc_v2/reads/DELETED
FILES=()
while read -r size name; do
  got=$(stat -c %s "$SRC/$name" 2>/dev/null || echo MISSING)
  [ "$got" = "$size" ] || { echo "SIZE MISMATCH $name: expected $size got $got"; exit 1; }
  FILES+=("$W/$name")
done < <(awk '$NF ~ /dc\.q20\.fastq$/ {print $5, $NF}' "$MANIFEST")
echo "== ${#FILES[@]} inputs match the syng.tsv set"

build() { # name, make args...
  local name=$1; shift
  mkdir "$W/src_$name"
  git -C "$SYNG" archive "$COMMIT" | tar -x -C "$W/src_$name"
  make -C "$W/src_$name" -j1 "$@" syng >/dev/null 2>&1   # -j1: PR5's syncmer_iter.o rule races
  cp "$W/src_$name/syng" "$W/syng_$name"
}
build seqhash AVX2=1
build csyncmer AVX2=1 CSYNCMER=1
echo "== built"

t0=$(date +%s)
for f in "${FILES[@]}"; do cp "$SRC/$(basename "$f")" "$f"; done
echo "== copied to tmpfs in $(( $(date +%s) - t0 ))s"; df -h /dev/shm | tail -1

# extraction symbols: seqhash's iterator, or csyncmer_fast + its syng wrappers
EXTRACT_RE='^(syncmerNext|syncmerIteratorReinit|syncmerIterator|seqhashIterator|syncmerMultiRead|csyncmer_.*)$'

for b in seqhash csyncmer; do
  rm -f "$W"/o.*
  /usr/bin/time -f "TIME %e user %U sys %S rss %M" taskset -c 0-7 \
    "$W/syng_$b" -w 992 -k 31 -T 8 -o "$W/o" "${FILES[@]}" > "$OUT/${b}.stdout" 2> "$OUT/${b}.time"
  echo "$b plain: $(grep ^TIME "$OUT/${b}.time")  $(grep -o 'instances of [0-9]* syncmers' "$OUT/${b}.stdout")"

  rm -f "$W"/o.*
  /usr/bin/time -f "TIME %e user %U sys %S rss %M" taskset -c 0-7 \
    perf record -F 499 -o "$W/perf_$b.data" \
    "$W/syng_$b" -w 992 -k 31 -T 8 -o "$W/o" "${FILES[@]}" > "$OUT/${b}_perf.stdout" 2> "$OUT/${b}_perf.time"
  echo "$b perf:  $(grep ^TIME "$OUT/${b}_perf.time")"
  perf report -i "$W/perf_$b.data" --no-children --sort sym --stdio --percent-limit 0 2>/dev/null \
    | grep -v '^#' | grep -v '^$' > "$OUT/perf_${b}_sym.txt" || true
  perf report -i "$W/perf_$b.data" --no-children --sort pid --stdio 2>/dev/null \
    | grep -v '^#' | grep -v '^$' > "$OUT/perf_${b}_thread.txt" || true
  cp "$W/perf_$b.data" "$OUT/"
  awk -v re="$EXTRACT_RE" '{ p=$1; sub(/%/,"",p); s=$NF; tot+=p; if (s ~ re) ex+=p }
       END { printf "extraction share %.2f%% of samples (sum of all listed %.1f%%)\n", ex, tot }' \
    "$OUT/perf_${b}_sym.txt" | tee "$OUT/${b}_share.txt"
  sed -n 1,15p "$OUT/perf_${b}_sym.txt"
  echo "-- per thread"; cat "$OUT/perf_${b}_thread.txt"
done
echo "== done $(date)"
