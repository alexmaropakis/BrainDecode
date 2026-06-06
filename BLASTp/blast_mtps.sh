#!/bin/bash
#SBATCH --job-name=blast_mtps
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/maropakis.a/logs/June26_blast_mtps_%j.out
#SBATCH --error=/scratch/maropakis.a/logs/June26_blast_mtps_%j.err

set -euo pipefail

export PATH=$HOME/bin/ncbi-blast-2.17.0+/bin:$PATH
blastp -version

THREADS=${SLURM_CPUS_PER_TASK:-16}
DEP=/scratch/maropakis.a/Dependencies

python3 /scratch/maropakis.a/scripts/BLAST_mtps.py \
  --csv        "$DEP/mtp_maps/SAAP_quant_df.csv" \
  --human-ref  "$DEP/FASTA/HUMAN.fasta" \
  --mouse-ref  "$DEP/FASTA/MOUSE_UP000000589_10090.fasta" \
  --out-dir    "$DEP/mtp_maps/" \
  --prefix     April26 \
  --threads    "$THREADS"
