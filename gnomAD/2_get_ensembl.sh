#!/bin/bash
# 2_get_ensembl.sh

# Fetch the Ensembl release-110 (GRCh38) protein FASTA and GTF, build the BLAST
# protein database, and regenerate gene_coords.tsv from the SAME GTF so gene
# symbols are release-consistent with the proteome used for BLAST in stage 4.

# Outputs (all under $ENSDIR):
#   Homo_sapiens.GRCh38.pep.all.fa            BLAST target (stage 4)
#   Homo_sapiens.GRCh38.pep.all.fa.p*         BLAST database index
#   Homo_sapiens.GRCh38.110.gtf               transcript / gene model
#   gene_coords.tsv                           gene -> chrom/start/end (stage 5)

#SBATCH --job-name=get_ensembl110
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/get_ensembl_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/get_ensembl_%j.err

set -euo pipefail

ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
REL=110
ENSDIR=$ROOT/ensembl_${REL}
mkdir -p "$ENSDIR" "$ROOT/logs"
cd "$ENSDIR"

BASE=https://ftp.ensembl.org/pub/release-${REL}/fasta/homo_sapiens
GTFBASE=https://ftp.ensembl.org/pub/release-${REL}/gtf/homo_sapiens

# 1. proteome (BLAST target). Headers carry transcript:ENST... gene_symbol:...
wget -c "$BASE/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

# 2. GTF (gene-region bounds for stage-5 scanning; source for gene_coords.tsv)
wget -c "$GTFBASE/Homo_sapiens.GRCh38.${REL}.gtf.gz"

gunzip -kf Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip -kf Homo_sapiens.GRCh38.${REL}.gtf.gz

# 3. gene_coords.tsv: gene_name -> chrom, start, end (1-based, from gene rows).
#    Regenerated from the Ensembl 110 GTF so symbols match the BLAST proteome.
#    Ensembl GTF uses bare contig names (1, 2, ..., X, Y); stage 5 maps these to
#    the gnomAD chr-prefixed names at query time.
gawk -F'\t' '$3=="gene"{
    if (match($9, /gene_name "([^"]+)"/, a))
        print a[1]"\t"$1"\t"$4"\t"$5
}' "Homo_sapiens.GRCh38.${REL}.gtf" > gene_coords.tsv

# 4. BLAST protein database from the Ensembl proteome
export PATH=$HOME/bin/ncbi-blast-2.17.0+/bin:$PATH
makeblastdb -in Homo_sapiens.GRCh38.pep.all.fa -dbtype prot -parse_seqids

echo "done: Ensembl ${REL} proteome, GTF, gene_coords.tsv, BLAST db in $ENSDIR"



#!/bin/bash
# 2_get_ensembl.sh 

# Fetch Ensembl GRCh38-110 release protein FASTA and GTF, 
# build BLAST protein database, 
# generate gene_coords.tsv from GTF

# OUTPUTS:
#   Homo_sapiens.GRCh38.pep.all.fa            BLAST target 
#   Homo_sapiens.GRCh38.pep.all.fa.p*         BLAST database index
#   Homo_sapiens.GRCh38.110.gtf               transcript / gene model
#   gene_coords.tsv                           gene -> chrom/start/end (stage 5)

#SBATCH --job-name=get_ensembl110
#SBATCH --partition=short
#SBATCH --time=24:00:00
#SBATCH --output=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/get_ensembl_%j.out
#SBATCH --error=/scratch/maropakis.a/Dependencies/gnomAD_pipeline/logs/get_ensembl_%j.err

set -euo pipefail

ROOT=/scratch/maropakis.a/Dependencies/gnomAD_pipeline
ENSDIR = $ROOT/ensembl_110
mkdir -p "$ENSDIR" "$ROOT/logs"
cd "$ENSDIR"

BASE=https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens
GTFBASE=https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens

# get proteome 
wget -c "$BASE/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

# get GTF (gene-region bounds for gene_coords.tsv)
wget -c "$GTFBASE/Homo_sapiens.GRCh38.${REL}.gtf.gz"

# Unpackage!!!
gunzip -kf Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip -kf Homo_sapiens.GRCh38.${REL}.gtf.gz

# Build gene_coords.tsv 
# maps gene_name to chromosome, start, end (1-based, from gene rows)
# generated from Ensembl 110 GTF, maps bare contig names (1,2,...,X,Y)
# t gnomAD chr-prefixed names at query time (basically different formatting needed)
gawk -F'\t' '$3=="gene"{
    if (match($9, /gene_name "([^"]+)"/, a))
        print a[1]"\t"$1"\t"$4"\t"$5
}' "Homo_sapiens.GRCh38.${REL}.gtf" > gene_coords.tsv

# build BLAST protein database from Ensembl proteome 
export PATH=$HOME/bin/ncbi-blast-2.17.0+/bin:$PATH
makeblastdb -in Homo_sapiens.GRCh38.pep.all.fa -dbtype prot -parse_seqids