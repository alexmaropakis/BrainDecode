#!/bin/bash
# 4_blast_bpseq.sh

# Run stage 4: BLAST base peptides against the Ensembl 110 proteome and write
# saap_mapped.csv (exact hits only, transcript IDs carried through).

#SBATCH --job-name=blast_mtps
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/mtp_maps/logs/BLASTp_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/mtp_maps/logs/BLASTp_%j.err

set -euo pipefail

ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
SCRIPTS=/home/maropakis.a/scripts/gnomAD
ENS=$ROOT/ensembl_110/Homo_sapiens.GRCh38.pep.all.fa

export PATH=$HOME/bin/ncbi-blast-2.17.0+/bin:$PATH
blastp -version

THREADS=${SLURM_CPUS_PER_TASK:-16}
mkdir -p "$ROOT/mtp_maps/logs"

python3 "$SCRIPTS/4_blast_bpseq.py" \
  --csv     "$ROOT/mtp_maps/saap_input.csv" \
  --ref     "$ENS" \
  --out     "$ROOT/mtp_maps/" \
  --threads "$THREADS"
