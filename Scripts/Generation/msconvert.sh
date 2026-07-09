#!/usr/bin/env bash
#SBATCH --job-name=raw2mzml
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/maropakis.a/spectra/logs/raw2mzml_%j.log

set -uo pipefail

TRFP=~/thermoRawFileParser/ThermoRawFileParser
PLEX=/scratch/maropakis.a/spectra/cortex_1_tsumagari # change to the plex you want to convert

echo "Processing plex: $PLEX"

for raw in "$PLEX"/*.raw; do
    [ -e "$raw" ] || continue
    base=$(basename "$raw" .raw)
    out="$PLEX/$base.mzML"

    # Skip only if output exists AND is a complete indexed mzML
    if [ -s "$out" ] && tail -c 200 "$out" | grep -q "</indexedmzML>"; then
        echo "SKIP $base"
        continue
    fi

    "$TRFP" -i="$raw" -o="$PLEX" -f=2 -l=3      # f=2 = indexed mzML

    if [ -s "$out" ] && tail -c 200 "$out" | grep -q "</indexedmzML>"; then
        echo "OK $base"
    else
        echo "FAIL $base"
    fi
done

# Submit:  sbatch 5_msconvert.sh
