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
ENSDIR=$ROOT/ensembl_110
mkdir -p "$ENSDIR" "$ROOT/logs"
cd "$ENSDIR"

BASE=https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens
GTFBASE=https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens

# get proteome 
wget -c "$BASE/pep/Homo_sapiens.GRCh38.pep.all.fa.gz"

# get GTF (gene-region bounds for gene_coords.tsv)
wget -c "$GTFBASE/Homo_sapiens.GRCh38.110.gtf.gz"

# Unpackage!!!
gunzip -kf Homo_sapiens.GRCh38.pep.all.fa.gz
gunzip -kf Homo_sapiens.GRCh38.110.gtf.gz

# Build gene_coords.tsv 
# maps gene_name to chromosome, start, end (1-based, from gene rows)
# generated from Ensembl 110 GTF, maps bare contig names (1,2,...,X,Y)
# t gnomAD chr-prefixed names at query time (basically different formatting needed)
gawk -F'\t' '$3=="gene"{
    if (match($9, /gene_name "([^"]+)"/, a))
        print a[1]"\t"$1"\t"$4"\t"$5
}' "Homo_sapiens.GRCh38.110.gtf" > gene_coords.tsv

# build BLAST protein database from Ensembl proteome 
export PATH=$HOME/bin/ncbi-blast-2.17.0+/bin:$PATH
makeblastdb -in Homo_sapiens.GRCh38.pep.all.fa -dbtype prot -parse_seqids