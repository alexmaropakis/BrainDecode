# gnomAD AF/SFS annotation pipeline

Annotate SAAP (substituted amino-acid peptide) sites with gnomAD v4.1 exome
allele frequencies and Ensembl 110 codon context, following the allele-frequency
analysis of Tsour et al. 2024 (Fig. 5j / Extended Data Fig. 10d). gnomAD-only and
Ensembl-only: no Roulette, no AlphaMissense, no constraint modelling.

## Layout

All scripts live in one folder:

    /home/maropakis.a/scripts/gnomAD/
        1_get_exome.sh
        2_get_ensembl.sh
        3_build_saap_csv.py   3_build_saap_csv.sh
        4_blast_bpseq.py      4_blast_bpseq.sh
        5_gnomad_parse.py     5_gnomad_parse.sh

All data lives under one root:

    ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
        ensembl_110/          proteome, GTF, gene_coords.tsv, BLAST db   (stage 2)
        mtp_maps/             saap_input.csv, saap_mapped.csv            (stages 3,4)
        gnomad/               <prefix>_gnomad_annotated.csv              (stage 5)
        logs/                 SLURM logs

    gnomAD VCFs (large, kept separate):
        /scratch/maropakis.a/gnomad_v4.1_exomes/

## Run order

    sbatch 1_get_exome.sh        # gnomAD v4.1 exome VCFs (skip if already present)
    sbatch 2_get_ensembl.sh      # Ensembl 110 proteome + GTF + gene_coords.tsv + BLAST db
    sbatch 3_build_saap_csv.sh   # pickles -> saap_input.csv   (set PREFIX inside)
    sbatch 4_blast_bpseq.sh      # saap_input.csv -> saap_mapped.csv
    sbatch 5_gnomad_parse.sh     # saap_mapped.csv -> <prefix>_gnomad_annotated.csv

Stages 1 and 2 are independent and can run in parallel. Stage 3 needs the AAS
pipeline pickles (validation2 + quantification already run). Stages 3->4->5 are
sequential.

## Data flow

    pickles --3--> saap_input.csv ((source_dataset, bp_seq, mtp_seq))
            --4--> saap_mapped.csv ((+ gene, transcript, protein_pos, bp_aa, mtp_aa))
            --5--> <prefix>_gnomad_annotated.csv

## Why Ensembl 110

Stage 4 BLASTs base peptides against the Ensembl 110 proteome (the release used
by Tsour et al. for substitution-context mapping) and keeps only 100%-identity,
zero-mismatch, zero-gap, full-length hits, so each protein position is exact and
the matched transcript (ENST) is known. gene_coords.tsv is regenerated from the
SAME Ensembl 110 GTF, so gene symbols are release-consistent across BLAST
(stage 4) and the region lookup (stage 5).

The on-disk HUMAN.fasta (UniProt) and HUMAN_GENOME.fna (NCBI RefSeq contigs)
are NOT used: UniProt headers carry no ENST anchor, and RefSeq contig names
(NC_000001.11) do not match an Ensembl GTF.

## Stage 5 output columns

Carried through from stage 4:
    dataset, gene, transcript, protein, protein_pos, bp_aa, mtp_aa, bp_seq, mtp_seq

Added from gnomAD / Ensembl CSQ:
    gnomad_variant       chr-pos-ref-alt of the matching variant (reproducing sub)
    gnomad_chr/pos/ref/alt
    gnomad_ac            allele count (AC)
    gnomad_af            allele frequency (AF); may be 0 (monomorphic)
    polymorphic          AF > 0 (only when the sub is reproduced by a variant)
    codon                reference codon (from VEP Codons)
    alt_codon            variant codon (only when the sub is reproduced)
    codon_change_pos     1-based position of the changed base in the codon
    ref_base / alt_base  reference / alternate nucleotide
    mutated_base         alternate nucleotide (alias of alt_base)
    strand               transcript strand (VEP STRAND)
    consequence          VEP Consequence
    reproduces_saap_sub  the gnomAD variant produces the SAME substitution
    exact_mutation       a single-nt change of the ref codon encodes mtp_aa
    gnomad_status        see below

### gnomad_status values

    reproduces_saap_sub      gnomAD has a variant making the SAME substitution;
                             AF/AC/polymorphic/variant populated
    site_variant_other_sub   gnomAD has missense at the site but a different
                             substitution; codon context populated, AF blank
                             (allele frequency is variant-specific)
    site_no_gnomad_missense  no observed missense at the site
    gene_not_found           gene symbol absent from gene_coords.tsv
    vcf_missing              expected per-chrom VCF not present
    fetch_error              tabix fetch failed (bad contig / index)
    unmappable_row          stage-4 row missing required fields

## Notes

- Annotate-everything, filter-never: no row is dropped; subset downstream on
  gnomad_af / polymorphic / exact_mutation / gnomad_status.
- AF=0 (monomorphic) is reported as 0 and is distinct from blank (no matching
  variant) -- the two are separable for downstream filtering.
- exact_mutation is computed only where gnomAD observed a missense (a CSQ codon
  exists to read). Sites gnomAD never observed get blank exact_mutation; to
  define it at every SAAP site, the codon must be read from the Ensembl CDS
  independently (not part of this AF/SFS path).
- Contig join: Ensembl gene coords are bare (1, 2, ..., X, Y); gnomAD VCF
  contigs are chr-prefixed (chr1). Stage 5 adds the prefix at query time.
