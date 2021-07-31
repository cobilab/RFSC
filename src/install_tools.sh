#!/bin/bash

Check_Installation() {
    if ! [ -x "$(command -v $1)" ]; then
        echo -e "\033[6;31mERROR\033[0;31m: $1 is not installed!" >&2;
        exit 1;
    else
        echo -e "\033[1;32mSUCCESS!\033[0m";
    fi
}

echo "Start Installation..."

conda install -c cobilab falcon --yes
conda install -c cobilab gto --yes
conda install -y -c bioconda geco3
conda install -c bioconda ac --yes
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
./src/install_orf.sh
#
pip install wheel
pip install pandas
pip install numpy
pip install sklearn
pip install keras
pip install tensorflow
#
Check_Installation "FALCON";
Check_Installation "gto";
Check_Installation "GeCo3";
Check_Installation "AC";
Check_Installation "cryfa";
Check_Installation "trimmomatic";
Check_Installation "fastp";
Check_Installation "spades.py";
Check_Installation "fastq_pair";
Check_Installation "efetch";
Check_Installation "blastn";
Check_Installation "dnadiff";
Check_Installation "art_illumina";
Check_Installation "./ORFs/ORFfinder";