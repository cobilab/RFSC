#!/bin/bash
#
testing_levels=(1 2 3 4 7 12 14);
#
VIRAL=0;
BACTERIA=0;
ARCHAEA=0;
PROTOZOA=0;
FUNGI=0;
PLANT=0;
MITOCHONDRION=0;
PLASTID=0;
HELP=0;
#
# ==============================================================================
#
if [ "$#" -eq 0 ];
  then
  HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help)                                          HELP=1; shift  ;;
    -vi|--viral)                VIRAL=1;                HELP=0; shift  ;;
    -ba|--bacteria)             BACTERIA=1;             HELP=0; shift  ;;
    -ar|--archaea)              ARCHAEA=1;              HELP=0; shift  ;;
    -pr|--protozoa)             PROTOZOA=1;             HELP=0; shift  ;;
    -fu|--fungi)                FUNGI=1;                HELP=0; shift  ;;
    -pl|--plant)                PLANT=1;                HELP=0; shift  ;;
    -mi|--mitochondrion )       MITOCHONDRION=1;        HELP=0; shift  ;;
    -ps|--plastid)              PLASTID=1;              HELP=0; shift  ;;
    *) echo "Invalid arg "$1 ; exit 1;
    ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
#
# ==============================================================================
# HELP
#
if [[ "$HELP" -eq "1" ]]
  then
  echo "                                                                                "
  echo "Usage: ./dna_compression.sh [options]                                           "
  echo "                                                                                "
  echo "It compresses and generates a CSV file with several levels of                   "
  echo "compression for the choosen NCBI database:                                      "
  echo "                                                                                "
  echo "   -h,   --help                   Show this help message,                       "
  echo "                                                                                "
  echo "   -vi,  --viral                  Genome compression viral DB,                  "
  echo "   -ba,  --bacteria               Genome compression bacteria DB,               "
  echo "   -ar,  --archaea                Genome compression archaea DB,                "
  echo "   -pr,  --protozoa               Genome compression protozoa DB,               "
  echo "   -fu,  --fungi                  Genome compression fungi DB,                  "
  echo "   -pl,  --plant                  Genome compression plant DB,                  "
  echo "   -mi,  --mitochondrion          Genome compression mitochondrion DB,          "
  echo "   -ps,  --plastid                Genome compression plastid DB.                "
  echo "                                                                                "
  echo "Example: ./dna_compression.sh --viral --fungi --plant                           "
  echo "                                                                                "
  exit 1
  fi
#
# ==============================================================================
#
################## VIRUS ##################
#
function VIRUS_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Virus Genomes"

    mkdir GeCo3_Output/Virus

    # for file in /media/alexloure/T7Touch/NCBI-Virus/NM-viral/*
    for file in References/NCBI-Virus/NM-viral/*
    do
        if [[ "$file" == "viral_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            printf "$out_file\n" >> GeCo3_Output/Virus/VIRAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Virus/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Virus/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/geco_VIRAL.txt
                rm GeCo3_Output/Virus/TO_COMPRESS*
            done
        fi
    done
}
#
################## BACTERIA ##################
#
function BACTERIA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Bacteria Genomes"

    mkdir GeCo3_Output/Bacteria

    for file in /media/alexloure/T7Touch/NCBI-Bacteria/NM-bacteria/*
    # for file in References/NCBI-Bacteria/NM-bacteria/*
    do
        if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            printf "$out_file\n" >> GeCo3_Output/Bacteria/BACTERIA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Bacteria/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Bacteria/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/geco_BACTERIA.txt
                rm GeCo3_Output/Bacteria/TO_COMPRESS*
            done
        fi
    done   
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Archaea Genomes"

    mkdir GeCo3_Output/Archaea

    # for file in /media/alexloure/T7Touch/NCBI-Archaea/NM-archaea/*
    for file in References/NCBI-Archaea/NM-archaea/*
    do
        if [[ "$file" == "archaea_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            printf "$out_file\n" >> GeCo3_Output/Archaea/ARCHAEA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Archaea/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Archaea/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/geco_ARCHAEA.txt
                rm GeCo3_Output/Archaea/TO_COMPRESS*
            done
        fi
    done
}
#
# ################## FUNGI ##################
#
function FUNGI_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Fungi Genomes"

    mkdir GeCo3_Output/Fungi

    # for file in /media/alexloure/T7Touch/NCBI-Fungi/NM-fungi/*
    for file in References/NCBI-Fungi/NM-fungi/*
    do
        if [[ "$file" == "fungi_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            printf "$out_file\n" >> GeCo3_Output/Fungi/FUNGI_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Fungi/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Fungi/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/geco_FUNGI.txt
                rm GeCo3_Output/Fungi/TO_COMPRESS*
            done
        fi
    done
}
#
# ################## PLANT ##################
#
function PLANT_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plant Genomes"

    mkdir GeCo3_Output/Plant

    # for file in /media/alexloure/T7Touch/NCBI-Plant/NM-plant/*
    for file in References/NCBI-Plant/NM-plant/*
    do
        if [[ "$file" == "plant_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            printf "$out_file\n" >> GeCo3_Output/Plant/PLANT_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plant/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Plant/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/geco_PLANT.txt
                rm GeCo3_Output/Plant/TO_COMPRESS*
            done
        fi
    done
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Protozoa Genomes"

    mkdir GeCo3_Output/Protozoa

    # for file in /media/alexloure/T7Touch/NCBI-Protozoa/NM-protozoa/*
    for file in References/NCBI-Protozoa/NM-protozoa/*
    do
        if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            printf "$out_file\n" >> GeCo3_Output/Protozoa/PROTOZOA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Protozoa/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Protozoa/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/geco_PROTOZOA.txt
                rm GeCo3_Output/Protozoa/TO_COMPRESS*
            done
        fi
    done
}
#
################## PLASTID ##################
#
function PLASTID_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start parsing Plastid Genomes"
    ulimit -n 1000000
    path="/media/alexloure/T7Touch/NCBI-Plastid"
    gunzip -k $path/DB-plastid.fa.gz
    mkdir $path/NM-plastid
    awk '/>/ {filename="/media/alexloure/T7Touch/NCBI-Plastid/NM-plastid/"NR".fna"} {print > filename}' $path/DB-plastid.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plastid Genomes"

    mkdir GeCo3_Output/Plastid

    for file in /media/alexloure/T7Touch/NCBI-Plastid/NM-plastid/*
    # for file in References/NCBI-Plastid/NM-plastid/*
    do
        if [[ "$file" == "plastid_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna"}"
            printf "$out_file\n" >> GeCo3_Output/Plastid/PLASTID_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plastid/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Plastid/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plastid/geco_PLASTID.txt
                rm GeCo3_Output/Plastid/TO_COMPRESS*
            done
        fi
    done
    #rm DB-plastid.*
}
#
# ################## MITOCHONDRIAL ##################
#
function MITOCHONDRIAL_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start parsing Mitochondrial Genomes"
    ulimit -n 1000000
    path="/media/alexloure/T7Touch/NCBI-Mitochondrial"
    gunzip -k $path/DB-mitochondrion.fa.gz
    mkdir $path/NM-mitochondrion
    awk '/>/ {filename="/media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/"NR".fna"} {print > filename}' $path/DB-mitochondrion.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

    mkdir GeCo3_Output/Mitochondrial

    # for file in References/NCBI-Mitochondrial/NM-mitochondrion/*
    for file in /media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/*
    do
        if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna"}"
            printf "$out_file\n" >> GeCo3_Output/Mitochondrial/MITOCHONDRIAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Mitochondrial/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Mitochondrial/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Mitochondrial/geco_MITOCHONDRIAL.txt
                rm GeCo3_Output/Mitochondrial/TO_COMPRESS*
            done
        fi
    done
    #rm DB-mitochondrion.*
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Viral DB"
  ./src/dna_compression_analysis_csv.sh Virus VIRAL geco3_Viral.csv
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_GENOME_COMPRESSION;
    echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Bacteria DB"
  ./src/dna_compression_analysis_csv.sh Bacteria BACTERIA geco3_Bacteria.csv
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Archeae DB"
  ./src/dna_compression_analysis_csv.sh Archaea ARCHAEA geco3_Archaea.csv
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Protozoa DB"
  ./src/dna_compression_analysis_csv.sh Protozoa PROTOZOA geco3_Protozoa.csv
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Fungi DB"
  ./src/dna_compression_analysis_csv.sh Fungi FUNGI geco3_Fungi.csv
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plant DB"
  ./src/dna_compression_analysis_csv.sh Plant PLANT geco3_Plant.csv
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Mitochondrial DB"
  ./src/dna_compression_analysis_csv.sh Mitochondrial MITOCHONDRIAL geco3_Mitochondrial.csv
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plastid DB"
  ./src/dna_compression_analysis_csv.sh Plastid PLASTID geco3_Plastid.csv
  fi
#
# ==============================================================================
#