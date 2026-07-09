#!/bin/bash
# 3_build_saap_csv.sh

# Run stage 3: extract MTP/BP pairs from the AAS pipeline pickles into
# saap_input.csv. Set PREFIX for the dataset you are processing (Ping_, ACG_,
# FC_, ...). Run once per prefix if you need separate inputs.

#SBATCH --job-name=build_saap
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/build_saap_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/build_saap_%j.err

set -euo pipefail

ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
SCRIPTS=/home/maropakis.a/scripts/gnomAD
AAS_ROOT=/scratch/maropakis.a/AAS_Pipeline
PREFIX=Ping_ # set per dataset

mkdir -p "$ROOT/mtp_maps/logs" "$ROOT/logs"

python3 "$SCRIPTS/3_build_saap_csv.py" \
  --root   "$AAS_ROOT/" \
  --prefix "$PREFIX" \
  --out    "$ROOT/mtp_maps/saap_input.csv"
