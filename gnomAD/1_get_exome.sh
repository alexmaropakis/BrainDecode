#!/bin/bash
# 1_get_exome.sh
#
# Download the gnomAD v4.1 exome site VCFs (one per chromosome) plus their
# tabix indices. These carry the allele frequencies (AF), allele counts (AC),
# and the embedded Ensembl VEP CSQ that the rest of the pipeline reads.
#
# Output: $DEST/gnomad.exomes.v4.1.sites.chr<N>.vcf.bgz (+ .tbi)

#SBATCH --job-name=get_gnomad
#SBATCH --partition=short
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/get_exome_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/get_exome_%j.err

set -euo pipefail

ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
DEST=/scratch/maropakis.a/gnomad_v4.1_exomes         
BASE=https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes

mkdir -p "$DEST" "$ROOT/logs"

CHROMS="1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y"
for c in $CHROMS; do
  f="gnomad.exomes.v4.1.sites.chr${c}.vcf.bgz"
  wget -c -P "$DEST" "$BASE/$f"
  wget -c -P "$DEST" "$BASE/$f.tbi"
done

echo "done: gnomAD v4.1 exomes in $DEST"
