#!/bin/bash
# 5_gnomad_parse.sh

# Run stage 5: annotate saap_mapped.csv against gnomAD v4.1 exomes via the
# embedded VEP CSQ. Output is the full annotated table; no filtering.

#SBATCH --job-name=gnomad_annot
#SBATCH --partition=short
#SBATCH --cpus-per-task=6
#SBATCH --mem=8G
#SBATCH --time=48:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/gnomad/logs/gnomad_annot_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/gnomad/logs/gnomad_annot_%j.err

set -euo pipefail

ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
SCRIPTS=/home/maropakis.a/scripts/gnomAD
REL=110

mkdir -p "$ROOT/gnomad/logs"

python3 "$SCRIPTS/5_gnomad_parse.py" \
  --keep-csv   "$ROOT/mtp_maps/saap_mapped.csv" \
  --gnomad-dir /scratch/maropakis.a/gnomad_v4.1_exomes \
  --gene-bed   "$ROOT/ensembl_${REL}/gene_coords.tsv" \
  --out-dir    "$ROOT/gnomad/" \
  --prefix     Ping # dataset name or whatever you want
