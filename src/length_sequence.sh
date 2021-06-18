#!/bin/bash

file=$1 #References/NCBI-Virus/GB_DB_VIRAL/GCA_000838505.1_ViralProj14028_genomic.fna.gz
dna_aa=$2
zipped=$3

if [[ "$dna_aa" == "DNA" ]]; then
    if [[ "$zipped" = "ZIPPED" ]]; then
        zcat $1 | grep -v ">" | tr -d -c "ACGT" > GENOME_FILE
    elif [[ "$zipped" = "UNZIPPED" ]]; then
        cat $1 > GENOME_FILE
    fi
elif [[ "$dna_aa" == "AA" ]]; then
    if [[ "$zipped" = "ZIPPED" ]]; then
        zcat $1 | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > GENOME_FILE
    elif [[ "$zipped" = "UNZIPPED" ]]; then
        cat $1 > GENOME_FILE
    fi
fi

SYMBOLS=`wc -c GENOME_FILE | awk '{print $1}'`

echo "$SYMBOLS"

rm GENOME_FILE