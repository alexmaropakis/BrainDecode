# gnomAD Annotation Pipeline

Annotate SAAP (substituted amino-acid peptide) substitution sites against
**gnomAD v4.1 exome** allele frequencies and **Ensembl 110** codon context.

The purpose is to ask, for every detected amino-acid substitution, **whether the
substitution that mass-spectrometry called as a mistranslation/misincorporation
event is also seen in the human population as a genetically encoded
polymorphism**. A SAAP site that coincides with a common missense SNP is
genomically explainable (it may be a real coding variant in the sample, not a
ribosomal error); a SAAP site that gnomAD never observed, or observed only as a
monomorphic/ultra-rare variant, is more consistent with a true translational
substitution. This pipeline produces the per-site evidence table needed to make
that distinction; it **does not** filter the dataset itself.

---

## 1. Inputs / Outputs

**Consumes:** the per-SAAP quantification pickles (`MTP_quant_dict.p`) produced by
the upstream AAS pipeline (validation2 + quantification must already have run),
i.e. the MTP/BP sequence pairs that define each substitution.

**Produces:** `<prefix>_gnomad_annotated.csv` — one row per substitution site,
carrying the original MS-derived identity plus gnomAD allele frequency/count,
Ensembl codon context, the substitution direction check, and a categorical
`gnomad_status` that records exactly how (or whether) gnomAD resolved the site.

**Scope:** gnomAD has no mouse data, so Takasugi (aged mouse) SAAPs
are not annotated here and pass downstream unannotated; this pipeline is run only on
the human datasets (Ping 2018 — ACG/FC plexes).

---

## 2. Repository layout

All scripts live in one folder:
```
/home/maropakis.a/scripts/gnomAD/
    1_get_exome.sh
    2_get_ensembl.sh
    3_build_saap_csv.py    3_build_saap_csv.sh
    4_blast_bpseq.py       4_blast_bpseq.sh
    5_gnomad_parse.py      5_gnomad_parse.sh
```
All derived data lives under one root:
```
ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
    ensembl_110/   proteome FASTA, GTF, gene_coords.tsv, BLAST db   (stage 2)
    mtp_maps/      saap_input.csv, saap_mapped.csv                  (stages 3, 4)
    gnomad/        <prefix>_gnomad_annotated.csv                    (stage 5)
    logs/          SLURM logs
```
gnomAD VCFs are large and kept separate:
```
/scratch/maropakis.a/gnomad_v4.1_exomes/
    gnomad.exomes.v4.1.sites.chr<N>.vcf.bgz  (+ .tbi)
```
All experimental data (MTP_quant_dict.p) is kept under a common root. 
{PREFIX} is the general name beginning every experimental folder, it's important for allowing the scripts to walk through every data directory to get the right .p files. Make sure to change {PREFIX} in each script and the paths to match your own. 
```
    AAS_ROOT=/scratch/maropakis.a/AAS_Pipeline
    PREFIX=Ping_ # set per dataset
        ├── Ping_2018_ACGb1 # example below of AAS Pipeline output files 
        │   ├── DP_dict.p
        │   ├── DP_search_evidence_dict.p
        │   ├── Ion_validated_MTP_dict.p
        │   ├── MTP_dict.p
        │   ├── ***MTP_quant_dict.p***
        │   ├── PTM_dict.p
        │   ├── qMTP_dict.p
        │   ├── tonsil_q-value_Precision_Recall_data.xlsx
        │   ├── Validated_MTP_dict.p
        │   └── Validation_search_evidence_dict.p
        ├── Ping_2018_ACGb2
        ├── Ping_2018_ACGb3
        ├── Ping_2018_ACGb4
        ├── Ping_2018_ACGb5
        ├── Ping_2018_FCb1
        ├── Ping_2018_FCb2
        ├── Ping_2018_FCb3
        ├── Ping_2018_FCb4
        ├── Ping_2018_FCb5
```
---

## 3. Pipeline

```
    sbatch 1_get_exome.sh        # gnomAD v4.1 exome VCFs
    sbatch 2_get_ensembl.sh      # Ensembl 110 proteome + GTF + gene_coords.tsv + BLAST db
    sbatch 3_build_saap_csv.sh   # pickles -> saap_input.csv   (set PREFIX inside)
    sbatch 4_blast_bpseq.sh      # saap_input.csv -> saap_mapped.csv
    sbatch 5_gnomad_parse.sh     # saap_mapped.csv -> <prefix>_gnomad_annotated.csv
```

Stages 1 and 2 are **independent** (downloads/reference build) and can run in
parallel. Stage 3 requires the AAS pipeline pickles. Stages **3 → 4 → 5 are
strictly sequential** (each consumes the previous output). Run stage 3 once per
dataset prefix if you want separate input tables.

---

## 4. Data flow

    pickle files (MTP_quant_dict.p)
      --3-->  saap_input.csv     (source_dataset, saap_id, bp_seq, mtp_seq, pickle_file)
      --4-->  saap_mapped.csv    (+ gene, transcript[ENST], protein, protein_pos, bp_aa, mtp_aa)
      --5-->  <prefix>_gnomad_annotated.csv   (+ gnomAD AF/AC, codon context, status)

Each stage adds the coordinate information the next stage needs: stage 4 turns
peptide sequences into **exact protein positions on a known transcript**, and
stage 5 uses that transcript+position to read gnomAD's embedded VEP annotation
directly, with no genomic-coordinate reconstruction required.

---

## 5. Stages of running the pipeline

### Stage 1 — `1_get_exome.sh` : fetch gnomAD v4.1 exome VCFs

Downloads the per-chromosome gnomAD v4.1 **exome** site VCFs (chr1–22, X, Y) and
their tabix indices from the GCP public bucket
(`storage.googleapis.com/gcp-public-data--gnomad/release/4.1/vcf/exomes`) into
`/scratch/maropakis.a/gnomad_v4.1_exomes/`. `wget -c` resumes partial downloads,
so a re-run is OK if it crashes.

These VCFs carry the allele frequency (`AF`), allele count (`AC`), and the
**embedded Ensembl VEP consequence annotation** that the rest of the pipeline
reads. The contigs are **`chr`-prefixed** (`chr1`, …), which matters for the
contig-name join in stage 5.

> *Caveat to verify before stage 5:* the VEP INFO key can be `vep` or `CSQ`
> depending on release. Stage 5 currently reads the `vep` INFO key. Confirm with
> `bcftools view -h <one VCF> | grep -E 'ID=(vep|CSQ)'`; if the exome release
> exposes `CSQ`, update the two `info["vep"]` / `info.get("vep")` references in
> `5_gnomad_parse.py`.

### Stage 2 — `2_get_ensembl.sh` : build the Ensembl 110 reference set

Downloads the Ensembl GRCh38 **release-110** protein FASTA
(`Homo_sapiens.GRCh38.pep.all.fa`) and GTF, then produces three artifacts:

1. **BLAST protein database** (`makeblastdb -dbtype prot -parse_seqids`) — the
   target for stage 4 base-peptide mapping.
2. **`gene_coords.tsv`** — `gene_name → (chrom, start, end)`, 1-based, parsed
   from `gene` rows of the GTF. Used by stage 5 to define the VCF scan region per
   gene. Contig names are **bare** (`1, 2, …, X, Y`), matching Ensembl space; the
   `chr` prefix is added at query time in stage 5.
3. The unzipped **proteome FASTA** kept on disk as the BLAST input.

Regenerating `gene_coords.tsv` from the *same* GTF that
backs the BLAST proteome keeps gene symbols release-consistent across stage 4 and
stage 5. The UniProt `HUMAN.fasta` and NCBI RefSeq `HUMAN_GENOME.fna` are
deliberately **not** used: UniProt headers carry no ENST anchor, and RefSeq
contig names (`NC_000001.11`) do not match an Ensembl GTF.

### Stage 3 — `3_build_saap_csv.py` / `.sh` : extract MTP/BP pairs from pickles

Walks experiment folders under `--root` whose names start with `--prefix`
(`Ping_`, `ACG_`, `FC_`, …), loads each `MTP_quant_dict.p`, and writes one CSV row
per SAAP entry. The TMT plex identifier (e.g. `ACGb1`, `FCB2`) is inferred from
the folder name via regex; entries missing `BP_seq` or `MTP_seq` are skipped.

**Output `saap_input.csv`:**

| column | meaning |
|---|---|
| `source_dataset` | TMT plex identifier inferred from the folder |
| `saap_id` | key index in the quant dictionary |
| `bp_seq` | base peptide sequence (the canonical/encoded peptide) |
| `mtp_seq` | mutant/substituted tryptic peptide sequence |
| `pickle_file` | origin file path (traceability) |

This is purely an **extraction/consolidation** step — no biology, just the
peptide pairs that become the BLAST queries.

### Stage 4 — `4_blast_bpseq.py` / `.sh` : map base peptides to exact protein coordinates

BLASTs each **base peptide** (`bp_seq`) against the Ensembl 110 proteome
(`blastp-short`, `-comp_based_stats 0`, `-evalue 1000`,
`-max_target_seqs 5`) and keeps **only fully exact hits**: 100% identity, 0
mismatches, 0 gaps, full peptide length, query anchored at residue 1
(`qstart==1`, `qend==qlen`). Among qualifying hits the best bitscore is taken.

Because the BP maps perfectly and full-length onto a protein, the protein
position of every substituted residue is **exact**:
`protein_pos = subject_start + query_offset`. The matched **transcript (ENST)** is
parsed from the Ensembl FASTA header (`transcript:` token) and carried through —
this is the key that lets stage 5 read gnomAD's per-transcript VEP annotation
directly. A header lacking a `transcript:` token is treated as a fatal format
mismatch (wrong proteome) and counted in a warning.

The script emits **one row per differing residue** between `bp_seq` and
`mtp_seq` (a SAAP can carry more than one substitution).

**Output `saap_mapped.csv`:**

| column | meaning |
|---|---|
| `dataset` | carried from `source_dataset` |
| `gene` | Ensembl gene symbol (from `gene_symbol:` token) |
| `transcript` | Ensembl transcript ID (ENST, versioned) |
| `protein` | Ensembl protein ID (ENSP) |
| `protein_pos` | 1-based protein position of the substituted residue |
| `bp_aa` | encoded amino acid (base-peptide residue) |
| `mtp_aa` | substituted/misincorporated amino acid |
| `bp_seq`, `mtp_seq` | the originating peptide pair |

Mapping uses base peptides (not the substituted peptides) because the BP is the
genome-encoded sequence and so must match a reference protein exactly; the
substitution is then read off as the position where BP and MTP differ. The
**direction is BP → MTP** throughout (encoded → misincorporated), consistent with
the rest of BrainDecode.

### Stage 5 — `5_gnomad_parse.py` / `.sh` : annotate against gnomAD

For each mapped site (keyed on **ENST + protein position + encoded aa**), scans
the gnomAD exome VCF over the gene's region (one tabix scan per gene; rows grouped
by gene for efficiency) and reads the embedded **Ensembl VEP CSQ**. Because
gnomAD's CSQ already supplies, per transcript, the `Codons`, the `Amino_acids`
change, the `Protein_position`, the `STRAND`, and the `Feature` (ENST), **no
genomic-coordinate reconstruction is needed** for any variant gnomAD observed.
Only `PASS`/unfiltered records and `missense_variant` consequences are considered;
`AC`/`AF` are read per matched alternate allele, keeping the max-AF record per key.

**Match logic (two-tier):**

1. **Preferred — substitution reproduced.** A gnomAD variant exists with, for the
   same alt allele, `Feature == ENST`, `Protein_position == protein_pos`,
   encoded aa `== bp_aa`, **and** incorporated aa `== mtp_aa`. This is the
   population variant making the *same* substitution as the SAAP. Its AF/AC and
   full variant identity are reported. → `gnomad_status = reproduces_saap_sub`.
2. **Fallback — site has other missense.** No variant reproduces the SAAP
   substitution, but gnomAD observed *some* missense at the site (encoded aa
   `== bp_aa`, different incorporated aa). Only the **codon context** is borrowed
   (codon, strand, exact-mutation check); **AF/AC are left blank** because allele
   frequency is variant-specific and must not be borrowed from a different
   substitution. → `gnomad_status = site_variant_other_sub`.

`exact_mutation` is an independent check: does a **single-nucleotide neighbour**
of the reference codon (which must encode `bp_aa`) translate to `mtp_aa`? It
flags whether the observed substitution is reachable by one SNV (BP → MTP
direction), and is only computable where a CSQ codon exists to read.

---

## 6. Stage 5 output columns (`<prefix>_gnomad_annotated.csv`)

**Carried through from stage 4:**

    dataset, gene, transcript, protein, protein_pos, bp_aa, mtp_aa, bp_seq, mtp_seq

**Added from gnomAD / Ensembl CSQ:**

| column | meaning |
|---|---|
| `gnomad_variant` | `chr-pos-ref-alt` of the matching variant (only when sub reproduced) |
| `gnomad_chr` / `gnomad_pos` / `gnomad_ref` / `gnomad_alt` | genomic coordinates of that variant |
| `gnomad_ac` | allele count (AC) |
| `gnomad_af` | allele frequency (AF); **may be 0** (monomorphic — meaningful, not missing) |
| `polymorphic` | `AF > 0` (set only when the sub is reproduced; blank otherwise) |
| `codon` | reference codon (from VEP `Codons`) |
| `alt_codon` | variant codon (only when sub reproduced) |
| `codon_change_pos` | 1-based position of the changed base within the codon |
| `ref_base` / `alt_base` | reference / alternate nucleotide |
| `mutated_base` | alternate nucleotide (alias of `alt_base`) |
| `strand` | transcript strand (VEP `STRAND`) |
| `consequence` | VEP `Consequence` |
| `reproduces_saap_sub` | the gnomAD variant produces the SAME substitution |
| `exact_mutation` | a single-nt change of the ref codon encodes `mtp_aa` |
| `gnomad_status` | categorical resolution status (below) |

### `gnomad_status` values

| value | meaning |
|---|---|
| `reproduces_saap_sub` | gnomAD has a variant making the **same** substitution; AF/AC/`polymorphic`/variant populated |
| `site_variant_other_sub` | gnomAD has missense at the site but a **different** substitution; codon context populated, AF blank (variant-specific) |
| `site_no_gnomad_missense` | no observed missense at this site in gnomAD |
| `gene_not_found` | gene symbol absent from `gene_coords.tsv` |
| `vcf_missing` | expected per-chromosome VCF not present |
| `fetch_error` | tabix fetch failed (bad contig / index) |
| `unmappable_row` | stage-4 row missing required fields (bad/short aa, non-numeric position, etc.) |

The status column is the heart of the design: it lets the analyst separate
**unresolvable coordinates** (`gene_not_found`, `vcf_missing`, `fetch_error`,
`unmappable_row`) from **genuinely resolved-but-absent** sites
(`site_no_gnomad_missense`) and from **resolved-and-present** sites
(`reproduces_saap_sub`, `site_variant_other_sub`). These are scientifically
different statements and are never collapsed.

---

## 7. Design principles

- **Annotate-everything, filter-never.** No input row is dropped; every row is
  emitted with a status. Subset downstream on `gnomad_af`, `polymorphic`,
  `exact_mutation`, or `gnomad_status`. This keeps pre-filter counts available and
  avoids baking analyst decisions into the reference table.
- **`AF = 0` is a distinct category, not "missing".** A monomorphic variant
  (observed in gnomAD's reference build but with zero observed alternate alleles)
  is reported as `0` and is separable from a blank AF (no matching variant). The
  two carry different meaning for the SFS / population-genetics interpretation.
- **AF is variant-specific and never borrowed.** When only a *different*
  substitution is observed at the site, codon context is reported but AF is left
  blank — attaching another variant's frequency would be wrong.
- **Direction is BP → MTP everywhere** (encoded → misincorporated), matching the
  corrected `compute_swap`/`swap` convention in the rest of BrainDecode.
- **Exact coordinates, exact matches only.** Stage 4 keeps only 100%/0mm/0gap/
  full-length BLAST hits, so positions are never approximate.
- **Contig-name join handled explicitly.** Ensembl gene coords are bare
  (`1, 2, …`); gnomAD VCF contigs are `chr`-prefixed (`chr1`). Stage 5 adds the
  prefix at query time, the single point where the two namespaces meet.
- **`exact_mutation` only where a codon exists.** It is computed only where gnomAD
  observed a missense (a CSQ codon is available). Sites gnomAD never observed get
  blank `exact_mutation`; defining it at *every* SAAP site would require reading
  the codon from the Ensembl CDS independently — outside this AF/SFS path.

---

## 8. Downstream utility

The annotated table is the genomic-plausibility layer for the SAAP set. Typical
uses:

1. **Separate encoded variants from translational substitutions.** Sites with
   `reproduces_saap_sub` and `polymorphic == True` (especially at appreciable
   `gnomad_af`) are candidates for being genuine coding polymorphisms in the
   sample rather than misincorporation events — i.e. a known confound to flag or
   remove before interpreting a SAAP as a ribosomal error. Sites that are
   `site_no_gnomad_missense`, or `reproduces_saap_sub` with `gnomad_af` at/near 0,
   are the population-absent / population-rare substitutions most consistent with
   true mistranslation.

2. **Site-frequency-spectrum analysis (Tsour-style).** The `gnomad_af` /
   `gnomad_ac` distribution across reproduced sites supports the AF/SFS panels
   that compare the substitution set against the population allele-frequency
   spectrum, distinguishing rare/ultra-rare from common variation.

3. **Single-SNV reachability.** `exact_mutation` indicates which substitutions are
   reachable by a single nucleotide change — relevant to whether a substitution is
   a plausible DNA-level event vs. a translational one, and to codon/near-cognate
   misreading interpretation.

4. **Joining back to quantification.** The table joins back to `SAAP_quant_df`
   so RAAS values can be analyzed conditional on gnomAD status — e.g. testing 
   whether aging-associated RAAS trajectories differ between population-absent 
   and population-present sites.

5. **Reporting categories.** `gnomad_status` provides ready-made buckets for
   supplemental tables (e.g. counts of reproduced vs. other-missense vs.
   absent-vs-unresolvable sites per dataset) without any additional bookkeeping.

---

## 9. Caveats

- **Human only.** gnomAD has no mouse; Takasugi mouse SAAPs are not annotated by
  this pipeline.
- **`max-target-seqs` interaction.** Stage 4 caps at 5 targets and then filters to
  exact full-length hits; if a base peptide is shared across many transcripts,
  confirm the intended transcript survives the exact-hit filter.
- **Codon context vs. AF.** For `site_variant_other_sub`, treat `codon`/`strand`/
  `exact_mutation` as valid but `gnomad_af`/`gnomad_ac` as intentionally blank.
- **Re-runnability.** Downloads (`wget -c`) and `makeblastdb` are idempotent; stage
  5 rewrites its output table each run.
