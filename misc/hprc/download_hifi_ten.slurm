#!/usr/bin/env bash
#SBATCH -J hprc10_dl
#SBATCH -p seqbio
#SBATCH -c 12
#SBATCH --mem=16G
#SBATCH -t 24:00:00
#SBATCH -o /pasteur/appa/scratch/rchikhi/hprc_v2/reads/logs/dl_%j.out
# Re-fetch the 10 HPRC Release 2 DeepConsensus files behind Fig 1(c) / syng.tsv
# and decompress them to the .fastq names syng_benchmark.sh expects.
# URLs resolved from human-pangenomics/hprc_intermediate_assembly
#   data_tables/sequencing_data/data_hifi_release2_v1.0.index.csv
set -uo pipefail
OUT=/pasteur/appa/scratch/rchikhi/hprc_v2/reads
S3=https://human-pangenomics.s3.amazonaws.com/submissions
A=42AFCE59-2972-45EE-BB72-5B84E7B57377--HPRC_DEEPCONSENSUS_v1pt2_2023_08_q20
B=3A25CF8A-1F36-42EE-BC9F-D29CECAA2A99--HPRC_DEEPCONSENSUS_v1pt2_2023_12_q20
FILES=(
  "$A HG00099 m54329U_220825_174247"  "$A HG00280 m54329U_220901_221341"
  "$A HG00558 m54329U_220107_233847"  "$A HG00639 m54329U_211222_104516"
  "$A HG01074 m54329U_211110_112322"  "$B HG01123 m54329U_200205_002609"
  "$B HG02257 m64076_200125_231256"   "$B HG02486 m64076_200211_192227"
  "$A HG02922 m54329U_220816_182601"  "$A HG03209 m64076_220526_115049"
)
fetch() {
  local sub=$1 s=$2 cell=$3
  local n="$s.$cell.dc.q20.fastq" url="$S3/$sub/$s/raw_data/PacBio_HiFi/deepconsensus/v1pt2/$s.$cell.dc.q20.fastq.gz"
  [ -s "$OUT/$n" ] && { echo "$n already present"; return 0; }
  for try in 1 2 3 4 5; do
    curl -sS --fail -L --retry 20 --retry-delay 10 -C - -o "$OUT/$n.gz.part" "$url" && break
    sleep 30
  done
  gzip -t "$OUT/$n.gz.part" || { echo "FAILED integrity: $n"; return 1; }
  pigz -dc -p 2 "$OUT/$n.gz.part" > "$OUT/$n.tmp" && mv "$OUT/$n.tmp" "$OUT/$n" && rm -f "$OUT/$n.gz.part"
  echo "done $n $(stat -c %s "$OUT/$n") bytes $(date)"
}
pids=()
for f in "${FILES[@]}"; do fetch $f & pids+=($!); done
rc=0; for p in "${pids[@]}"; do wait $p || rc=1; done
ls -la "$OUT"/*.fastq; du -cb "$OUT"/*.fastq | tail -1
exit $rc
