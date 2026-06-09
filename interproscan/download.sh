#!/bin/bash
#SBATCH --job-name=interproscan-download
#SBATCH --partition=short

mkdir InterProScan
cd InterProScan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.77-108.0/interproscan-5.77-108.0-64-bit.tar.gz
tar -pxvzf interproscan-5.77-108.0-*-bit.tar.gz

# where:
#     p = preserve the file permissions
#     x = extract files from an archive
#     v = verbosely list the files processed
#     z = filter the archive through gzip
#     f = use archive file
