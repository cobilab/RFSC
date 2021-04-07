#!/bin/bash

Check_Installation() {
    if ! [ -x "$(command -v $1)" ]; then
        echo -e "\e[41mERROR\e[49m: $1 is not installed!" >&2;
        exit 1;
    else
        echo -e "\e[42mSUCCESS!\e[49m";
    fi
}

echo "Start Installation..."

conda install -c cobilab falcon --yes
conda install -c cobilab gto --yes
conda install -y -c bioconda geco2
conda install -c bioconda -y cryfa
conda install -c bioconda trimmomatic --yes
conda install -c bioconda fastp --yes
conda install -c bioconda spades --yes
conda install -c bioconda fastq-pair --yes
conda install -c bioconda entrez-direct --yes
conda install -c bioconda/label/cf201901 entrez-direct --yes
conda install -c bioconda blast --yes
conda install -c bioconda mummer4 --yes
conda install -c bioconda art --yes
#
Check_Installation "FALCON";
Check_Installation "gto";
Check_Installation "GeCo2";
Check_Installation "cryfa";
Check_Installation "trimmomatic";
Check_Installation "fastp";
Check_Installation "spades.py";
Check_Installation "fastq_pair";
Check_Installation "efetch";
Check_Installation "blastn";
Check_Installation "dnadiff";
Check_Installation "art_illumina";