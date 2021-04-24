#!/bin/bash

file=$1 #References/NCBI-Virus/GB_DB_VIRAL/GCA_000838505.1_ViralProj14028_genomic.fna.gz

zcat $1 | grep -v ">" | tr -d -c "ACGT" > GENOME_FILE
grep -o 'GC' GENOME_FILE | wc -l
rm GENOME_FILE