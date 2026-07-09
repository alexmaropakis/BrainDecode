#!/usr/bin/env bash
# run_all_plexes.sh
#
# The ONE place that lists every TMT plex across all datasets. Each line is a self-contained
# call to gen_fragpipe_plex.py that prepares that plex end-to-end (annotation, raw->mzML,
# spectra staging, workflow+manifest) and writes submit/submit_<plex>.sh.
#
# ADDING A DATASET = add one block below (one gen line per plex). No other script edits anywhere.
# Channel count is auto-derived from each sample_map, so no -t needed.
#
# Disambiguated tokens (because tissue names repeat across studies):
#   Takasugi kidney -> 'kidney'      Keele kidney   -> 'kidney_keele'
#   Keele cortex    -> 'cortex_keele'  Tsumagari cortex rep1 -> 'cortex_1_tsumagari'
# These match the names already in annotations/ and sample_map/ on disk, and must match the
# plex token in each *_MTP.fasta / *_fragpipe.fasta (named S?_<token>_MTP.fasta).
#
# Prereqs per plex: a *_fragpipe.fasta built by 1_PrepFASTA.py + 2_buildFragFASTA.py, and a
# sample_map .xlsx. Plexes whose FASTA isn't built yet will error at the FASTA lookup -- build
# those first (see NOTE at bottom).
#
# NOT INCLUDED (not TMT plexes -- different pipeline):
#   PD_2026         : DIA multi-protease (GluC/LysC/Trypsin), DIA-Umpire SE -> MSFragger mass-offset
#   Giansanti_2022  : DIA (Giansanti_2022_DIA)
#
# TODO before running the new datasets: replace the TODO_*.workflow names below with your real
# templates (Bai / Keele / Tsumagari acquisition-specific workflows).

#SBATCH --job-name=buildplexes
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=32G
#SBATCH --time=48:00:00

set -euo pipefail

GEN=/home/maropakis.a/scripts/search_gen/Alex_gen_fragpipe.py
RAW=/scratch/maropakis.a/MQ_raw
SMAP=/scratch/maropakis.a/Dependencies/sample_map
FASTA=/scratch/maropakis.a/Dependencies/FASTA_fragpipe
TPL=/home/maropakis.a/scripts/search_gen/FragPipe/templates
OUT=/scratch/maropakis.a/Frag_outputs
SPECTRA=/scratch/maropakis.a/spectra

# gen <plex_token> <species> <workflow_file> <raw_subdir> <sample_map_file>
gen() {
  python3 "$GEN" "$RAW/$4" \
    --plex "$1" --species "$2" \
    --workflow "$TPL/$3" --sample-map "$SMAP/$5" \
    --fasta-dir "$FASTA" --out-dir "$OUT" --spectra-root "$SPECTRA"
}

# ===== Ping_2018  (human, TMT10 MS3): ACG b1-5, FC b1-5 =====
gen acgb1 human TMT10_MS3_Val.workflow Ping_2018/ACG/b1 sample_map_acgb1.xlsx
gen acgb2 human TMT10_MS3_Val.workflow Ping_2018/ACG/b2 sample_map_acgb2.xlsx
gen acgb3 human TMT10_MS3_Val.workflow Ping_2018/ACG/b3 sample_map_acgb3.xlsx
gen acgb4 human TMT10_MS3_Val.workflow Ping_2018/ACG/b4 sample_map_acgb4.xlsx
gen acgb5 human TMT10_MS3_Val.workflow Ping_2018/ACG/b5 sample_map_acgb5.xlsx
gen fcb1  human TMT10_MS3_Val.workflow Ping_2018/FC/b1  sample_map_fcb1.xlsx
gen fcb2  human TMT10_MS3_Val.workflow Ping_2018/FC/b2  sample_map_fcb2.xlsx
gen fcb3  human TMT10_MS3_Val.workflow Ping_2018/FC/b3  sample_map_fcb3.xlsx
gen fcb4  human TMT10_MS3_Val.workflow Ping_2018/FC/b4  sample_map_fcb4.xlsx
gen fcb5  human TMT10_MS3_Val.workflow Ping_2018/FC/b5  sample_map_fcb5.xlsx

# ===== Bai_2020  (human, pooled,TMT10 MS2) =====
gen pooled human TMT10_MS2_Val.workflow Bai_2020 sample_map_pooled.xlsx

# ===== Takasugi_2024  (mouse, TMT16 MS2): 8 tissues =====
gen aorta  mouse TMT16_Val.workflow Takasugi_2024/aorta  sample_map_aorta.xlsx
gen brain  mouse TMT16_Val.workflow Takasugi_2024/brain  sample_map_brain.xlsx
gen heart  mouse TMT16_Val.workflow Takasugi_2024/heart  sample_map_heart.xlsx
gen kidney mouse TMT16_Val.workflow Takasugi_2024/kidney sample_map_kidney.xlsx
gen liver  mouse TMT16_Val.workflow Takasugi_2024/liver  sample_map_liver.xlsx
gen lung   mouse TMT16_Val.workflow Takasugi_2024/lung   sample_map_lung.xlsx
gen muscle mouse TMT16_Val.workflow Takasugi_2024/muscle sample_map_muscle.xlsx
gen skin   mouse TMT16_Val.workflow Takasugi_2024/skin   sample_map_skin.xlsx

# ===== Keele_2025  (mouse, Astral TMT): TODO your workflow. tokens = <tissue>_keele =====
# gen cortex_keele      mouse TODO_KEELE.workflow Keele_2025/cortex      sample_map_cortex_keele.xlsx
# gen hippocampus_keele mouse TODO_KEELE.workflow Keele_2025/hippocampus sample_map_hippocampus_keele.xlsx
# gen kidney_keele      mouse TODO_KEELE.workflow Keele_2025/kidney      sample_map_kidney_keele.xlsx
# gen striatum_keele    mouse TODO_KEELE.workflow Keele_2025/striatum    sample_map_striatum_keele.xlsx

# ===== Tsumagari_2023  (mouse, TMT10 MS2): tokens = <tissue>_<rep>_tsumagari =====
# gen cortex_1_tsumagari      mouse TODO_TSUMAGARI.workflow Tsumagari_2023/cortex_1      sample_map_cortex_1_tsumagari.xlsx
gen cortex_2_tsumagari      mouse TMT10_MS2_Val.workflow Tsumagari_2023/cortex_2      sample_map_cortex_2_tsumagari.xlsx
gen hippocampus_1_tsumagari mouse TMT10_MS2_Val.workflow Tsumagari_2023/hippocampus_1 sample_map_hippocampus_1_tsumagari.xlsx
gen hippocampus_2_tsumagari mouse TMT10_MS2_Val.workflow Tsumagari_2023/hippocampus_2 sample_map_hippocampus_2_tsumagari.xlsx

echo
echo "All plexes prepped. Submit them with:"
echo "  for s in $OUT/submit/submit_*.sh; do sbatch \"\$s\"; done"

# NOTE: new-dataset FASTAs not built yet. Before the Bai/Keele/Tsumagari lines succeed you must
# run stages 1-2 for them: put their *_MTP.fasta in FASTA_appended/ named S?_<token>_MTP.fasta
# (e.g. S9_cortex_keele_MTP.fasta), pass their MQ roots to 1_PrepFASTA.py under --mouse-root /
# --human-root, then run 2_buildFragFASTA.py. Ping + Takasugi FASTAs already exist.
