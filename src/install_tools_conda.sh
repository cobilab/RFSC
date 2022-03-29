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
#
conda install -c cobilab falcon --yes
conda install -c cobilab gto --yes
conda install -y -c bioconda geco3
conda install -c bioconda ac --yes
conda install -c bioconda -y cryfa
conda install -c bioconda trimmomatic --yes
conda install -c bioconda fastp --yes
conda install -c bioconda fastq-pair --yes
conda install -c bioconda entrez-direct --yes
conda install -c bioconda/label/cf201901 entrez-direct --yes
conda install -c bioconda mummer4 --yes
conda install -c bioconda art --yes
conda install -c bioconda seqtk --yes
conda install -c conda-forge ncbi-datasets-cli --yes
conda install -c bioconda perl-lwp-protocol-https --yes
conda install -c bioconda sra-tools entrez-direct --yes
conda install perl-io-socket-ssl perl-net-ssleay perl-lwp-protocol-https entrez-direct --yes
conda install -c bioconda orfm --yes
#
./src/install_orf.sh
#
