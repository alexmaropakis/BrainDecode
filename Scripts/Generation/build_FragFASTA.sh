#!/bin/bash
#SBATCH --job-name=buildfragfasta
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time=02:00:00

# Run when you're ready to build the FragPipe FASTA files

set -euo pipefail

cd /home/maropakis.a/scripts/search_gen/ # run from the dir holding the .py scripts
DEP=/scratch/maropakis.a/Dependencies
mkdir -p "$DEP/mtp_maps" "$DEP/FASTA_fragpipe"

# --- Stage 1: Generate per-species CSVs (query id = full header + TMT tag) ---
  python3 1_PrepFASTA.py \
    --mtp-dir     /scratch/maropakis.a/Dependencies/FASTA_appended/ \
    --human-root  /scratch/maropakis.a/MQ_outputs/Ping_2018 \
    --human-root  /scratch/maropakis.a/MQ_outputs/Bai_2020 \
    --mouse-root  /scratch/maropakis.a/MQ_outputs/Takasugi_2024 \
    --mouse-root /scratch/maropakis.a/MQ_outputs/Tsumagari_2023 \
    --out-dir     /scratch/maropakis.a/Dependencies/mtp_maps/

# --- Stage 2: build per-plex FragPipe FASTAs from the per-plex CSVs ---
# each *_MTP.fasta uses ONLY its own {token}.csv; no cross-plex keep-set sharing
python3 2_buildFragFASTA.py \
  --mtp-dir "$DEP/FASTA_appended/" \
  --csv-dir "$DEP/mtp_maps/" \
  --out-dir "$DEP/FASTA_fragpipe/"

# --- Stage 3: sanity check — no duplicate headers in any FragPipe FASTA ---
fail=0
for f in "$DEP"/FASTA_fragpipe/*_fragpipe.fasta; do
  dups=$(grep '^>' "$f" | sort | uniq -d || true)
  if [ -n "$dups" ]; then
    echo "ERROR: duplicate headers in $f:" >&2
    echo "$dups" | head >&2
    fail=1
  fi
done
[ "$fail" -eq 0 ] && echo "All FragPipe FASTAs have unique headers." || exit 1
