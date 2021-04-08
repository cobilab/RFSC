#!/bin/bash

cd References/NCBI-Virus/

echo -e "\033[1;34m[RFSC]\033[0m Building the database at References/NCBI-Virus/";

wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > ASCG.txt
mkdir -p GB_DB_VIRAL
mkdir -p GB_DB_VIRAL_CDS
mkdir -p GB_DB_VIRAL_RNA
cat ASCG.txt | xargs -I{} -n1 -P8 wget -P GB_DB_VIRAL {}/*_genomic.fna.gz
mv GB_DB_VIRAL/*_cds_from_genomic.fna.gz GB_DB_VIRAL_CDS/
mv GB_DB_VIRAL/*_rna_from_genomic.fna.gz GB_DB_VIRAL_RNA/
zcat GB_DB_VIRAL/*.fna.gz > VDB.fa

echo -e "\033[1;34m[RFSC] \033[1;32m Building has been successful! \033[0m";