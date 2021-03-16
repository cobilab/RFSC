<p align="center"> <img src="Logo.png"> </p>

## 1. About ##

## 2. Reference databases ##

### Building a reference database from NCBI:

```
    wget ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/viral/assembly_summary.txt
    awk -F '\t' '{if($12=="Complete Genome") print $20}' assembly_summary.txt > ASCG.txt
    mkdir -p GB_DB_VIRAL
    mkdir -p GB_DB_VIRAL_CDS
    mkdir -p GB_DB_VIRAL_RNA
    cat ASCG.txt | xargs -I{} -n1 -P8 wget -P GB_DB_VIRAL {}/*_genomic.fna.gz
    mv GB_DB_VIRAL/*_cds_from_genomic.fna.gz GB_DB_VIRAL_CDS/
    mv GB_DB_VIRAL/*_rna_from_genomic.fna.gz GB_DB_VIRAL_RNA/
    zcat GB_DB_VIRAL/*.fna.gz > VDB.fa
```
