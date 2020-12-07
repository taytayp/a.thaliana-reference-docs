#!/usr/bin/env bash

# Taylor Falk
# tfalk@bu.edu
# processing RNASeq data for abundance levels

# create folder for storing everything
mkdir arabidopsis
cd arabidopsis

## step 1: trimmomatic
# download the paird end reads
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/002/ERR3333412/ERR3333412_1.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR333/002/ERR3333412/ERR3333412_2.fastq.gz

module load trimmomatic
# run trimmomatic using the compressed reads and Illumina-PE tags to filter out
trimmomatic PE -threads 6 ERR3333412_1.fastq.gz ERR3333412_2.fastq.gz out_forward_paired.fq.gz out_forward_unpaired.fq.gz out_reverse_paired.fq.gz out_reverse_unpaired.fq.gz ILLUMINACLIP:Illumina-PE.fasta:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36

## step 2: kallisto
# download Araport coding sequence reference fasta
wget https://www.arabidopsis.org/download_files/Genes/Araport11_genome_release/Araport11_blastsets/Araport11_genes.201606.cdna.fasta.gz
gzip -d Araport11_genes.201606.cdna.fasta.gz

# create an index with Araport genes
kallisto index -i Araport11_cds Araport11_genes.201606.cdna.fasta

mkdir kallisto 

# run kallisto to quantify and annotate rna seq data with gene names from Araport
kallisto quant -i Araport11_cds -o kallisto -t 6 out_forward_paired.fq.gz out_reverse_paired.fq.gz

# capture the protein names only
awk '{print $1}' abundance.tsv > rosette_leaf1_proteins.txt

# retrieve only the lines that do not have an ending zero, indicating abundance > 0
grep -P -v "0$" abundance.tsv > non-zero_abundance.tsv

# save these nonzero protein names
awk '{print $1}' non-zero_abundance.tsv > rosette_leaf1_proteins_nonzero.txt
