#!/bin/bash

# Set the genome version
genome_version="hg38"

# Create a directory to store the reference files
mkdir -p "${genome_version}_reference"
cd "${genome_version}_reference"

# Download the chromosome sizes file from UCSC
wget "http://hgdownload.soe.ucsc.edu/goldenPath/${genome_version}/bigZips/${genome_version}.chrom.sizes" -O "${genome_version}.chrom.sizes"

# Create a genome sizes file
awk -v OFS='\t' '{print $1, $2}' "${genome_version}.chrom.sizes" > "${genome_version}.genomeSizes_subset"

# Download the reference genome FASTA file from UCSC
wget "http://hgdownload.soe.ucsc.edu/goldenPath/${genome_version}/bigZips/${genome_version}.fa.gz"
gunzip "${genome_version}.fa.gz"

# Index the reference genome FASTA file using samtools
samtools faidx "${genome_version}.fa"

# Split the reference genome FASTA file into individual chromosome files
mkdir chr_fasta
awk 'BEGIN {RS = ">" ; FS = "\n" ; ORS = ""} {if ($2) {print ">"$0 > "chr_fasta/"$1".fa"}}' "${genome_version}.fa"

# Generate mappability tracks using gem-mappability
# You would need to install GEM (Genome Multi-tool) for this step
# Example command:
gem-mappability -I "${genome_version}.fa" -l 50 -o "${genome_version}_mappability"

# Navigate back to the original directory
cd ..
