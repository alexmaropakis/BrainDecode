#!/usr/bin/env bash
# run_all_plexes_DIA.sh
#
# The ONE place that lists every label-free DIA run across all datasets. Each line is a
# self-contained call to gen_fragpipe_labelfree_DIA.py that prepares that run end-to-end
# (spectra collection + workflow+manifest from metadata.csv) and writes submit/submit_<dataset>.sh.
#
# ADDING A DATASET = add one block below (one gen line per protease). No other script edits
# anywhere. Experiment/bioreplicate are auto-derived from each dataset's metadata.csv, so no
# channel/-t knobs -- every file becomes its own quant column grouped by Condition.
#
# This is the DIA sibling of run_all_plexes.sh (TMT). The differences that matter:
#   - no sample_map / channels; grouping comes from metadata.csv (ID,Condition,...).
#   - spectra are Bruker diaPASEF (.d/.d.dia) read DIRECTLY from their raw folders -- no
#     conversion, no staging into spectra/. FragPipe/diaTracer handles the .d natively.
#   - a dataset is split into one protease per subdir (PD_2026/{Trypsin,GluC,LysC}), and each
#     protease has its OWN .workflow template that already carries that protease's digest. So we
#     do NOT pass --enzyme -- the generator only patches db-path and leaves the digest alone.
#
# Dataset tokens must match the plex token in each *_fragpipe.fasta (named S?_<token>_fragpipe.fasta):
#   PD_2026 Trypsin -> 'pd_trypsin'   GluC -> 'pd_gluc'   LysC -> 'pd_lysc'
# (FASTAs are already built + named to this convention -- S1_pd_trypsin_fragpipe.fasta, etc.)
#
# Per-protease workflow templates live in $TPL (moved out of PD_2026/Trypsin/):
#   Trypsin -> DIA.workflow   GluC -> DIA_GluC.workflow   LysC -> DIA_LysC.workflow

#SBATCH --job-name=buildDIA
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail

GEN=/home/maropakis.a/scripts/search_gen/gen_fragpipe_labelfree_DIA.py
RAW=/scratch/maropakis.a/MQ_raw
FASTA=/scratch/maropakis.a/Dependencies/FASTA_fragpipe
TPL=/home/maropakis.a/scripts/search_gen/FragPipe/templates
OUT=/scratch/maropakis.a/Frag_outputs

# gen <dataset_token> <species> <workflow_file> <raw_subdir> <metadata_file>
#   <workflow_file> is relative to TPL and already carries that protease's digest (no --enzyme).
#   <metadata_file> is relative to RAW (usually the dataset's own metadata.csv, shared across
#   its proteases).
gen() {
  python3 "$GEN" "$RAW/$4" \
    --dataset "$1" --species "$2" \
    --workflow "$TPL/$3" --metadata "$RAW/$5" \
    --fasta-dir "$FASTA" --out-dir "$OUT"
}

# ===== PD_2026  (human diaPASEF, multi-protease; diaTracer -> MSFragger -> DIA-NN) =====
# One shared metadata.csv (PD_2026/metadata.csv); three proteases, three FASTAs, three runs.
gen pd_trypsin human DIA.workflow      PD_2026/Trypsin PD_2026/metadata.csv
gen pd_gluc    human DIA_GluC.workflow PD_2026/GluC    PD_2026/metadata.csv
gen pd_lysc    human DIA_LysC.workflow PD_2026/LysC    PD_2026/metadata.csv

# ===== Giansanti_2022  (human DIA, single protease): TODO your workflow =====
# gen giansanti_2022_dia human TODO_GIANSANTI.workflow Giansanti_2022 Giansanti_2022/metadata.csv

echo
echo "All label-free DIA runs prepped. Submit them with:"
echo "  for s in $OUT/submit/submit_*.sh; do sbatch \"\$s\"; done"
