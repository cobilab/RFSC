#!/bin/bash
# =============================================
# Instalation of other programs (without conda)
# =============================================
# Install ORFfinder
#
mkdir ORFs
cd ORFs
wget ftp://ftp.ncbi.nlm.nih.gov/genomes/TOOLS/ORFfinder/linux-i64/ORFfinder.gz
gunzip ORFfinder.gz 
chmod a+x ORFfinder 
wget https://github.com/wwood/OrfM/releases/download/v0.7.1/orfm-0.7.1_Linux_x86_64.tar.gz
tar xzf orfm-0.7.1_Linux_x86_64.tar.gz
mv orfm-0.7.1_Linux_x86_64 orfm
rm orfm-0.7.1_Linux_x86_64.tar.gz
chmod a+x orfm/orfm
#