#!/bin/bash

file=$1 #References/NCBI-Virus/GB_DB_VIRAL/GCA_000838505.1_ViralProj14028_genomic.fna.gz
zipped=$2

A=0;
T=0;
G=0;
C=0;

if [[ "$zipped" = "ZIPPED" ]]; then
    zcat $1 | grep -v ">" | tr -d -c "ACGT" > GENOME_FILE
elif [[ "$zipped" = "UNZIPPED" ]]; then
    cat $1 > GENOME_FILE
fi

A=`grep -o 'A' GENOME_FILE | wc -l`
T=`grep -o 'T' GENOME_FILE | wc -l`
G=`grep -o 'G' GENOME_FILE | wc -l`
C=`grep -o 'C' GENOME_FILE | wc -l`

GC=$(( $G + $C ))
ATGC=$(( $A + $T + $G + $C ))
PERCENTAGE=`echo "scale=8 ; $GC / $ATGC" | bc`

echo "GC: $GC ATGC: $ATGC PERCENTAGE: $PERCENTAGE"

rm GENOME_FILE