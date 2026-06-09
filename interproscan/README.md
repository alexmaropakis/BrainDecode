# BLASTp — MTP/SAAP Mapping & FragPipe FASTA Construction

Part of the [FragPipe Linux Tutorial](https://github.com/alexmaropakis/FragPipe_Linux_Tutorial). This step maps mistranslated peptides (MTPs / single–amino-acid substitution peptides) back to their parent proteins via `blastp`, classifies them, and builds FragPipe-ready FASTA databases that fold the validated MTP entries in alongside the reference proteome.

## Overview

The pipeline runs in three stages, orchestrated by a single SLURM script:

1. **BLAST mapping** (`blast_mtps.py`) — splits `*_MTP.fasta` inputs by species (human vs. mouse), blasts every MTP against the matching reference proteome, classifies each peptide, and writes per-species `_unfiltered.csv` and `_filtered.csv`.
2. **FASTA construction** (`build_fragpipe_fasta.py`) — uses the filtered CSVs to build FragPipe/Philosopher-compatible FASTAs: reference entries pass through unchanged, kept MTPs get mock-UniProt headers, and reversed decoys are appended for all targets.
3. **Sanity check** — verifies no duplicate headers exist in any output FASTA.

Species routing is filename-based throughout: tokens beginning with `ACG` or `FC` are treated as **human**; everything else (Aorta, Brain, Heart, …) is **mouse**.

## Requirements

- Linux environment with Bash and Python 3
- NCBI BLAST+ (`makeblastdb`, `blastp`) — the script was built against `ncbi-blast-2.17.0+`
- SLURM (for the provided batch script; the Python scripts can also be run standalone)
- Reference proteome FASTAs for human and mouse

## Input layout

- `FASTA_appended/` — input `*_MTP.fasta` files, one per TMT set (e.g. `ACG_B5_MTP.fasta`). MTP records are flagged by a `>MTP|` header prefix; non-MTP records are treated as reference entries.
- `FASTA/HUMAN.fasta`, `FASTA/MOUSE_UP000000589_10090.fasta` — reference proteomes used as BLAST databases.

## Usage

### Full pipeline (SLURM)

Submit from the directory holding the `.py` scripts:

```bash
sbatch build_fragfasta.sh
```

Edit the `#SBATCH` directives and the `DEP` path at the top of the script to match your environment. The script stages outputs into `Dependencies/mtp_maps/` and `Dependencies/FASTA_fragpipe/`.

### Running stages individually

**Stage 1 — BLAST and classify:**

```bash
python3 blast_mtps.py \
  --mtp-dir    Dependencies/FASTA_appended/ \
  --human-ref  Dependencies/FASTA/HUMAN.fasta \
  --mouse-ref  Dependencies/FASTA/MOUSE_UP000000589_10090.fasta \
  --out-dir    Dependencies/mtp_maps/ \
  --threads    16
```

**Stage 2 — Build FragPipe FASTAs:**

```bash
python3 build_fragpipe_fasta.py \
  --mtp-dir   Dependencies/FASTA_appended/ \
  --human-csv Dependencies/mtp_maps/human_filtered.csv \
  --mouse-csv Dependencies/mtp_maps/mouse_filtered.csv \
  --out-dir   Dependencies/FASTA_fragpipe/
```

## How it works

### BLAST parameters

Each MTP is searched with `blastp-short` (tuned for short query sequences), `-evalue 1000`, `-comp_based_stats 0`, and `-max_target_seqs 50`. Composition-based statistics are disabled and the e-value is deliberately permissive so that single-substitution matches are not filtered out before classification. Output is tabular (`-outfmt 6`) carrying `qseqid sseqid pident length mismatch gapopen qstart qend qlen stitle bitscore`. The reference DB is built on demand with `makeblastdb -parse_seqids`.

### Query IDs

Each query ID is the **full MTP header plus a TMT-set tag** (e.g. `...|TMT=acgb5`), with whitespace collapsed. This guarantees that the same MTP number reused across different TMT sets never collides in the query FASTA or in per-query hit grouping. Exact duplicates (same ID + same sequence) are skipped; same ID with a different sequence is uniquified with a `|dup2` suffix.

### Classification

A hit only counts as a full-length alignment if it has no gaps, starts at query position 1, and spans the entire query length. From those:

| Status | Meaning |
|---|---|
| `drop_reference` | An exact (0-mismatch) full-length match exists → peptide is just the reference, discard. |
| `keep` | Best full-length **1-mismatch** hit → genuine substitution; record accession + gene. |
| `drop_no_parent` | No 1-mismatch full-length parent found. |
| `drop_trypsin` | Best parent is trypsin (`PRSS1/2/3`, `TRY*`, or "trypsin" in title). |
| `drop_ig` | Best parent is immunoglobulin (`IGH/K/L`, `JCHAIN`, `IGJ`, or "immunoglobulin" in title). |

A `keep` entry is flagged `ambiguous` when its 1-mismatch candidate hits map to more than one gene.

### FragPipe FASTA output

For each input file, `build_fragpipe_fasta.py` writes `<name>_fragpipe.fasta` where:

- Reference entries pass through unchanged.
- Kept MTPs (matched by sequence against the filtered CSV) receive a Philosopher-safe mock-UniProt header with accession `<acc>-<mtp_id>-<tmt_set>`, made unique per entry.
- Reversed `rev_` decoys are appended for every target **except** the MTP entries (headers containing `-mut ` are skipped, so substitution peptides do not get decoys).

## Outputs

- `mtp_maps/{human,mouse}_unfiltered.csv` — every MTP entry with its classification.
- `mtp_maps/{human,mouse}_filtered.csv` — `keep` entries only.
- `mtp_maps/mtp_query_{species}.fasta`, `mtp_blast_{species}.tsv` — intermediate query and raw BLAST tables.
- `FASTA_fragpipe/*_fragpipe.fasta` — FragPipe-ready databases (reference + kept MTPs + decoys).

CSV columns: `query_id, mtp_id, tmt_set, source_file, sequence, species, status, accession, gene, ambiguous`.
