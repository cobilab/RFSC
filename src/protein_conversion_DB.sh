#!/bin/bash
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

ORFFINDER=0;
ORFM=1;
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
    -orffinder)                 ORFFINDER=1; ORFM=0;    HELP=0; shift  ;;
    -orfm)                      ORFM=1;                 HELP=0; shift  ;;
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
  echo "Usage: ./protein_conversion_DB.sh [options]                                     "
  echo "                                                                                "
  echo "Converts nucleotide sequences into protein sequences                            "
  echo "                                                                                "
  echo "   -h,   --help                   Show this help message,                       "
  echo "                                                                                "
  echo "   -orffinder                     Converion using ORFfinder                     "
  echo "   -orfm                          Convert using orfM (default)                  "
  echo "                                                                                "
  echo "   -vi,  --viral                  Convert viral DB,                             "
  echo "   -ba,  --bacteria               Convert bacteria DB,                          "
  echo "   -ar,  --archaea                Convert archaea DB,                           "
  echo "   -pr,  --protozoa               Convert protozoa DB,                          "
  echo "   -fu,  --fungi                  Convert fungi DB,                             "
  echo "   -pl,  --plant                  Convert plant DB,                             "
  echo "   -mi,  --mitochondrion          Convert mitochondrion DB,                     "
  echo "   -ps,  --plastid                Convert plastid DB.                           "
  echo "                                                                                "
  echo "Example: ./protein_conversion_DB.sh --viral --fungi --plant                     "
  echo "                                                                                "
  exit 1
  fi
#
# ==============================================================================
#
################## VIRUS ##################
#
function VIRUS_PROTEIN_COVERSION () {

    echo -e "\033[1;34m[RFSC]\033[0m Start Virus Protein Synthesizing"
    mkdir References/NCBI-Virus/PT-Virus

    for file in References/NCBI-Virus/NM-viral/*
    do
        if [[ "$file" == "viral_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"

            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > References/NCBI-Virus/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in References/NCBI-Virus/TO_PROTEIN -outfmt 0 -out References/NCBI-Virus/PROTEIN.fasta
                cat References/NCBI-Virus/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > References/NCBI-Virus/PT-Virus/$out_file.fna
                rm References/NCBI-Virus/TO_PROTEIN References/NCBI-Virus/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > References/NCBI-Virus/PROTEIN.fasta
                cat References/NCBI-Virus/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > References/NCBI-Virus/PT-Virus/$out_file.fna
                rm References/NCBI-Virus/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done
}
#
################## BACTERIA ##################
#
function BACTERIA_PROTEIN_COVERSION () {
    
    echo -e "\033[1;34m[RFSC]\033[0m Start Bacteria Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Bacteria/PT-Bacteria

    for file in /media/alexloure/T7Touch/NCBI-Bacteria/NM-bacteria/*
    # for file in References/NCBI-Bacteria/NM-bacteria/*
    do
        if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"

            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Bacteria/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Bacteria/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Bacteria/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Bacteria/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Bacteria/PT-Bacteria/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Bacteria/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Bacteria/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Bacteria/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Bacteria/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Bacteria/PT-Bacteria/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Bacteria/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done   
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_PROTEIN_COVERSION () {
    
    #already_converted=399;
    #stop_condition=1;

    echo -e "\033[1;34m[RFSC]\033[0m Start Archaea Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Archaea/PT-Archaea


    for file in /media/alexloure/T7Touch/NCBI-Archaea/NM-archaea/*
    #for file in References/NCBI-Archaea/NM-archaea/*
    do
        if [[ "$file" == "archaea_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            
            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Archaea/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Archaea/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Archaea/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Archaea/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Archaea/PT-Archaea/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Archaea/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Archaea/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Archaea/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Archaea/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Archaea/PT-Archaea/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Archaea/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done

}
#
# ################## FUNGI ##################
#
function FUNGI_PROTEIN_COVERSION () {
    
    echo -e "\033[1;34m[RFSC]\033[0m Start Fungi Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Fungi/PT-Fungi

    for file in /media/alexloure/T7Touch/NCBI-Fungi/NM-fungi/*
    #for file in References/NCBI-Fungi/NM-fungi/*
    do
        if [[ "$file" == "fungi_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            
            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Fungi/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Fungi/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Fungi/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Fungi/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Fungi/PT-Fungi/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Fungi/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Fungi/PROTEIN.fasta            
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Fungi/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Fungi/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Fungi/PT-Fungi/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Fungi/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done
}
#
# ################## PLANT ##################
#
function PLANT_PROTEIN_COVERSION () {
    
    echo -e "\033[1;34m[RFSC]\033[0m Start Plant Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Plant/PT-Plant

    for file in /media/alexloure/T7Touch/NCBI-Plant/NM-plant/*
    #for file in References/NCBI-Plant/NM-plant/*
    do
        if [[ "$file" == "plant_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"

            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Plant/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Plant/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Plant/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Plant/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Plant/PT-Plant/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Plant/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Plant/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Plant/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Plant/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Plant/PT-Plant/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Plant/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_PROTEIN_COVERSION () {
    
    echo -e "\033[1;34m[RFSC]\033[0m Start Protozoa Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Protozoa/PT-Protozoa

    for file in /media/alexloure/T7Touch/NCBI-Protozoa/NM-protozoa/*
    #for file in References/NCBI-Protozoa/NM-protozoa/*
    do
        if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna.gz"}"
            
            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Protozoa/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Protozoa/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Protozoa/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Protozoa/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Protozoa/PT-Protozoa/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Protozoa/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Protozoa/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Protozoa/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Protozoa/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Protozoa/PT-Protozoa/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Protozoa/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done
}
#
################## PLASTID ##################
#
function PLASTID_PROTEIN_COVERSION () {

    echo -e "\033[1;34m[RFSC]\033[0m Start Plastid Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Plastid/PT-Plastid

    for file in /media/alexloure/T7Touch/NCBI-Plastid/NM-plastid/*
    # for file in References/NCBI-Plastid/NM-plastid/*
    do
        if [[ "$file" == "plastid_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna"}"
            
            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Plastid/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Plastid/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Plastid/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Plastid/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Plastid/PT-Plastid/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Plastid/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Plastid/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Plastid/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Plastid/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Plastid/PT-Plastid/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Plastid/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done
}
#
# ################## MITOCHONDRIAL ##################
#
function MITOCHONDRIAL_PROTEIN_COVERSION () {

    echo -e "\033[1;34m[RFSC]\033[0m Start Mitochondrial Protein Synthesizing"
    mkdir /media/alexloure/T7Touch/NCBI-Mitochondrial/PT-Mitochondrial

    # for file in References/NCBI-Mitochondrial/NM-mitochondrion/*
    for file in /media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/*
    do
        if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        else
            out_file=$(basename $file);
            out_file="${out_file%".fna"}"
            
            if [[ "$ORFFINDER" -eq "1" ]]; then
                zcat $file | grep -v ">" | tr -d -c "ACGT" > /media/alexloure/T7Touch/NCBI-Mitochondrial/TO_PROTEIN
                ./ORFs/ORFfinder -s 1 -in /media/alexloure/T7Touch/NCBI-Mitochondrial/TO_PROTEIN -outfmt 0 -out /media/alexloure/T7Touch/NCBI-Mitochondrial/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Mitochondrial/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Mitochondrial/PT-Mitochondrial/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Mitochondrial/TO_PROTEIN /media/alexloure/T7Touch/NCBI-Mitochondrial/PROTEIN.fasta
            else
                ./ORFs/orfm/orfm $file > /media/alexloure/T7Touch/NCBI-Mitochondrial/PROTEIN.fasta
                cat /media/alexloure/T7Touch/NCBI-Mitochondrial/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > /media/alexloure/T7Touch/NCBI-Mitochondrial/PT-Mitochondrial/$out_file.fna
                rm /media/alexloure/T7Touch/NCBI-Mitochondrial/PROTEIN.fasta
            fi
            echo "Processed $out_file.fna"
        fi
    done
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_PROTEIN_COVERSION;
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_PROTEIN_COVERSION;
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_PROTEIN_COVERSION;
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_PROTEIN_COVERSION;
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_PROTEIN_COVERSION;
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_PROTEIN_COVERSION;
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_PROTEIN_COVERSION;
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_PROTEIN_COVERSION;
  fi
#
# ==============================================================================
#