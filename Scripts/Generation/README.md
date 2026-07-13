# `Scripts/Generation` — search & pipeline input generators

This folder holds the **generators** that turn raw mass-spec data + sample metadata into
ready-to-run search jobs (MaxQuant / FragPipe) and into the per-experiment scripts for the
SAAP detection pipeline (Tsour et al., *Nature* 2026). Nothing here does the analysis itself —
each script *writes* the config, FASTA, manifest, or SLURM job that a downstream tool then runs.

The design principle throughout: **adding a dataset should require no code edits** — you pass its
parameters on the command line (or add one line to a `run_all_*.sh` driver). Every script is
idempotent where it can be (skips existing mzML, existing frame-translation pickles, etc.).

---

## The big picture

```
                        ┌─────────────────────────────────────────────┐
  one-time reference    │  Alex_gen_translations.py                   │
                        │  genome → 6-frame translations + suffix arr │  (homology removal
                        └─────────────────────────────────────────────┘   for AAS pipeline)

  raw spectra ─┬─► MaxQuant route ───────────────────────────────────────────────────────┐
               │     Alex_gen_mqXML.py            (DP  = discovery / dependent-peptide)    │
               │     Alex_gen_Validation_mqXML.py (Val = validation vs *_MTP.fasta, v2.8)  │
               │     Alex_gen_v2880_Validation_mqXML.py (Val, MaxQuant v2.8.8.0)           │
               │                                                                           │
               └─► FragPipe route                                                          │
                     TMT (DDA):   Alex_gen_fragpipe.py           (one plex per call)       │
                                  Alex_gen_fragpipe_experiments.py (many plexes, one run)  │
                                  run_all_plexes.sh              (driver: every plex)      │
                     Label-free:  gen_fragpipe_labelfree.py      (DIA, raw→mzML)           │
                     (DIA)        gen_fragpipe_labelfree_DIA.py  (Bruker .d, in place)     │
                                  run_all_plexes_DIA.sh          (driver: every run)       │
                                                                                           │
       spectra conversion helpers:  Alex_msconvert.py / msconvert.sh  (raw → mzML)         │
                                                                                           ▼
  FragPipe FASTA build (MTP → FragPipe target+decoy FASTA):                                │
       1_PrepFASTA.py  →  2_buildFragFASTA.py   (driver: build_FragFASTA.sh)  ─────────────┘

  downstream SAAP pipeline step scripts (Detection → Validation1 → Validation2 → Quant):
       Alex_gen_pipeline.py   (fills path placeholders in step templates + writes SLURM)
```

`DP` = MaxQuant *dependent-peptide* search (first pass; finds the substituted/mistranslated
peptides). `Val` = validation search of those peptides (appended `*_MTP.fasta`) against the data.
`MTP` = mistranslated peptide; `SAAP` = substituted amino-acid peptide.

---

## Scripts by role

### Reference build (run once per genome)
| Script | What it produces |
|---|---|
| `Alex_gen_translations.py` | 6-frame genome translations + suffix-array pickles (`W{f}_aa_ambig.p`, `s{f}a_ambig.p`) for human & mouse, used by the AAS pipeline to strip homologous sequences. Paths are hard-coded in the `GENOMES` dict at the top; edit there to add a genome. Skips frames already on disk. |

### MaxQuant search generation
All three read a template mqpar `.xml`, rewrite the FASTA block / raw-file list / experiment &
fraction blocks / output folder, then write both the patched XML (`~/scripts/XML/<type>/`) and a
matching SLURM script (`~/scripts/Batch/searches/`).

| Script | Search | Notes |
|---|---|---|
| `Alex_gen_mqXML.py` | **DP** (discovery) | FASTA = species default (`HUMAN.fasta` / `MOUSE_…fasta`) or `--fasta-path`. Default MaxQuant `1.6.17.0`. Writes `MQ_outputs/<exp>/DP/<search_name>_DP/`. |
| `Alex_gen_Validation_mqXML.py` | **Val** | FASTA = the matching appended `*_MTP.fasta` (auto-matched from `--outfile` key for Ping/Takasugi). MaxQuant `2.8.0.0`, `MaxQuantCmd.exe`. |
| `Alex_gen_v2880_Validation_mqXML.py` | **Val** | Same idea, tuned for MaxQuant **v2.8.8.0** (`MaxQuantCmd.dll`, `OX=` taxonomy rule, stricter Takasugi FASTA suffix match). Use this when running the newer MQ build. |

### FragPipe FASTA build (MTP peptides → searchable FASTA)
Turns each appended `*_MTP.fasta` into a FragPipe-ready target+decoy FASTA, resolving every MTP to
its parent protein **strictly within its own plex** (no cross-plex leakage).

| Script | Step |
|---|---|
| `1_PrepFASTA.py` | For each plex, match every MTP peptide to its base peptide (single-residue diff) in that plex's MaxQuant **DP** `evidence.txt`, resolve the parent accession/gene/description from `proteinGroups.txt`, and write a per-plex `{token}.csv`. |
| `2_buildFragFASTA.py` | From each `{token}.csv` + its `*_MTP.fasta`, emit `S#_<token>_fragpipe.fasta`: references pass through, kept MTPs get a Philosopher-safe mock-UniProt header, and `rev_` decoys are appended for every target. Species tag is read from the CSV, never guessed. |
| `build_FragFASTA.sh` | SLURM driver: runs stage 1 → stage 2 → a duplicate-header sanity check over all output FASTAs. |

### FragPipe search generation — TMT (DDA)
| Script | Scope |
|---|---|
| `Alex_gen_fragpipe.py` | **One plex per call.** Builds `annotation.txt` from a single-plex sample map, converts `.raw`→`.mzML`, stages into `<spectra-root>/<plex>/`, patches the workflow (`db-path`, `tmtintegrator.channel_num`), and writes the manifest + `submit_<plex>.sh`. |
| `Alex_gen_fragpipe_experiments.py` | **One run, many experiments (plexes) → one manifest.** Each `--experiment NAME RAW_DIR PLEX` is a plex; `PLEX` selects its rows from a multi-plex sample map (matched against the `TMT plex` or `sample_ID` column). Stages each experiment into `<spectra-root>/<run>/<experiment>/` with its own `annotation.txt`, folds them all into one workflow + manifest + `submit_<run>.sh`. Point `--fasta` straight at the FASTA. |
| `run_all_plexes.sh` | The single place listing every TMT plex across all datasets — one `gen` line per plex calling `Alex_gen_fragpipe.py`. Add a dataset = add lines here. |

**Per-plex vs experiment-level:** use `Alex_gen_fragpipe.py` (via `run_all_plexes.sh`) when each
plex is its own independent search; use `Alex_gen_fragpipe_experiments.py` when several plexes must
be searched and quantified **together** in one FragPipe run (TMT-Integrator groups by experiment
across plexes).

### FragPipe search generation — label-free (DIA)
| Script | Scope |
|---|---|
| `gen_fragpipe_labelfree.py` | One label-free **DIA** dataset per call. No TMT/channels. Manifest built from a `metadata.csv` (`experiment=Condition`, per-file `bioreplicate`, `datatype=DIA`). Converts `.raw`→`.mzML` and stages; can patch the MSFragger enzyme per protease with `--enzyme`. |
| `gen_fragpipe_labelfree_DIA.py` | Same, but for **Bruker diaPASEF** (`.d`/`.d.dia`): spectra are read **in place** (no conversion, no staging), and each protease's workflow already carries its own digest (so `--enzyme` is an optional override only). |
| `run_all_plexes_DIA.sh` | Driver listing every label-free DIA run (one `gen` line per protease), calling `gen_fragpipe_labelfree_DIA.py`. |

### Spectra conversion helpers
| Script | Use |
|---|---|
| `Alex_msconvert.py` | Standalone `.raw`→`.mzML` (ThermoRawFileParser `-f=2`) + symlink into `<spectra-root>/<plex>/`. A thin subset of what the FragPipe TMT generators do inline — handy when you only need conversion/staging. |
| `msconvert.sh` | SLURM version for one plex directory, with an `</indexedmzML>` completeness check so partial conversions are re-run. Edit the `PLEX=` line to the folder you want. |

### Downstream SAAP pipeline step generation
| Script | Use |
|---|---|
| `Alex_gen_pipeline.py` | Generates the experiment-specific **Detection → Validation1 → Validation2 → Quant** scripts from templates in `~/scripts/Pipeline/templates`, filling in path placeholders (MQ dir, FragPipe results dir, sample map, FASTA base name) and writing one SLURM job per step. Filters templates by `--mode` (tmt/labelfree), `--species` (human/mouse), and `--engine` (maxquant/fragpipe). |

---

## Shared conventions

**Dataset / plex tokens** — lowercased, non-alphanumerics stripped for compact names, underscores
kept for disambiguated ones (tissue names repeat across studies):
`acgb1`, `fcb3`, `pooled`, `aorta`, `cortex_keele`, `cortex_1_tsumagari`, `pd_trypsin`.
The same token ties together a plex's sample map, its FASTA, and its outputs — keep it consistent.

**FASTA naming** (`Dependencies/FASTA_appended/` and `Dependencies/FASTA_fragpipe/`):
- `S#_<token>_noMTP.fasta` — reference (no substituted peptides)
- `S#_<token>_MTP.fasta` — reference **+** appended mistranslated peptides
- `S#_<token>_fragpipe.fasta` — FragPipe target+decoy build (from `2_buildFragFASTA.py`)

**On-disk layout** (`/scratch/maropakis.a/`):
```
MQ_raw/<Dataset>/<...>/            raw spectra (.raw, or Bruker .d/.d.dia)
MQ_outputs/<Dataset>/{DP,Val}/<search_name>_{DP,Val}/combined/txt/   MaxQuant results
Dependencies/FASTA/                reference FASTAs + genomes
Dependencies/FASTA_appended/       *_MTP.fasta / *_noMTP.fasta
Dependencies/FASTA_fragpipe/       *_fragpipe.fasta
Dependencies/mtp_maps/             per-plex {token}.csv from 1_PrepFASTA.py
Dependencies/sample_map/           per-plex TMT sample maps (.xlsx) + metadata.csv
Dependencies/frame_translations/   suffix-array pickles per species
spectra/<plex_or_run>/...          staged mzML + annotation.txt for FragPipe
Frag_outputs/{workflows,manifests,annotations,submit,logs,results}/  FragPipe run artifacts
```

**Sample maps** — TMT: an `.xlsx` with `TMT plex`, `TMT channel`, `sample_name`, `sample_ID`
columns (one plex, or several stacked for the multi-experiment generator). Label-free: a
`metadata.csv` with `ID` + `Condition` columns.

---

## Typical end-to-end sequences

**TMT via FragPipe**
1. MaxQuant DP search — `Alex_gen_mqXML.py` → run → `…/DP/…/combined/txt/`.
2. (AAS pipeline Detection produces the appended `*_MTP.fasta`.)
3. Build FragPipe FASTAs — `build_FragFASTA.sh` (`1_PrepFASTA.py` → `2_buildFragFASTA.py`).
4. Generate runs — `run_all_plexes.sh` (per plex) **or** `Alex_gen_fragpipe_experiments.py`
   (plexes that belong in one run), then `sbatch Frag_outputs/submit/submit_*.sh`.
5. Downstream — `Alex_gen_pipeline.py --engine fragpipe …` for Validation2/Quant.

**TMT via MaxQuant**
1. `Alex_gen_mqXML.py` (DP) → run.
2. `Alex_gen_Validation_mqXML.py` (or `…_v2880_…` for MQ 2.8.8.0) → run Val vs `*_MTP.fasta`.
3. `Alex_gen_pipeline.py --engine maxquant …`.

**Label-free DIA via FragPipe**
1. Build FASTAs (per protease) — `build_FragFASTA.sh`.
2. `run_all_plexes_DIA.sh` (Bruker `.d`) → `sbatch` the submit scripts.
3. `Alex_gen_pipeline.py --mode labelfree …`.

---

## Notes

- The two `gen_fragpipe_labelfree*.py` scripts and the `*.sh` drivers predate the `Alex_` naming
  convention used by the rest of the folder; their filenames are referenced verbatim inside the
  `run_all_*.sh` drivers, so renaming them means updating those `GEN=` paths too.
- `Alex_msconvert.py`'s `main()` references a bare `plex`/uses `a.raw_dir` loosely — it works for
  the common case but is the least-polished helper here; prefer the conversion built into the
  FragPipe TMT generators when you're staging for a search.
