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
#