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
  echo "Usage: ./protein_compression.sh [options]                                       "
  echo "                                                                                "
  echo "It compresses and generates a CSV file with several levels of                   "
  echo "compression for the choosen protein database:                                   "
  echo "                                                                                "
  echo "   -h,   --help                   Show this help message,                       "
  echo "                                                                                "
  echo "   -vi,  --viral                  Protein compression viral DB,                 "
  echo "   -ba,  --bacteria               Protein compression bacteria DB,              "
  echo "   -ar,  --archaea                Protein compression archaea DB,               "
  echo "   -pr,  --protozoa               Protein compression protozoa DB,              "
  echo "   -fu,  --fungi                  Protein compression fungi DB,                 "
  echo "   -pl,  --plant                  Protein compression plant DB,                 "
  echo "   -mi,  --mitochondrion          Protein compression mitochondrion DB,         "
  echo "   -ps,  --plastid                Protein compression plastid DB.               "
  echo "                                                                                "
  echo "Example: ./protein_compression.sh --viral --fungi --plant                       "
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

    if [[ ! -d "GeCo3_Output/Virus" ]]; then
        mkdir GeCo3_Output/Virus
    fi

    if [[ ! -d "GeCo3_Output/Virus/PT" ]]; then
        mkdir GeCo3_Output/Virus/PT
    fi

    if [[ -f "GeCo3_Output/Virus/PT/VIRAL_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Virus/PT/VIRAL_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Virus/PT-Virus/*
    for file in References/NCBI-Virus/PT-Virus/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "viral_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            
            printf "$out_file\n" >> GeCo3_Output/Virus/PT/VIRAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Virus/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Virus/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/PT/geco_VIRAL.txt
                rm GeCo3_Output/Virus/PT/TO_COMPRESS*
            done
        fi
    done
}
#
################## BACTERIA ##################
#
function BACTERIA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Bacteria Genomes"

    if [[ ! -d "GeCo3_Output/Bacteria" ]]; then
        mkdir GeCo3_Output/Bacteria
    fi

    if [[ ! -d "GeCo3_Output/Bacteria/PT" ]]; then
        mkdir GeCo3_Output/Bacteria/PT
    fi

    if [[ -f "GeCo3_Output/Bacteria/PT/BACTERIA_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Bacteria/PT/BACTERIA_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Bacteria/PT-Bacteria/*
    # for file in References/NCBI-Bacteria/PT-Bacteria/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Bacteria/PT/BACTERIA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Bacteria/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Bacteria/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/PT/geco_BACTERIA.txt
                rm GeCo3_Output/Bacteria/PT/TO_COMPRESS*
            done
        fi
    done   
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Archaea Genomes"

    if [[ ! -d "GeCo3_Output/Archaea" ]]; then
        mkdir GeCo3_Output/Archaea
    fi

    if [[ ! -d "GeCo3_Output/Archaea/PT" ]]; then
        mkdir GeCo3_Output/Archaea/PT
    fi

    if [[ -f "GeCo3_Output/Archaea/PT/ARCHAEA_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Archaea/PT/ARCHAEA_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Archaea/PT-Archaea/*
    #for file in References/NCBI-Archaea/PT-Archaea/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "archaea_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        
        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";
        
        else
            printf "$out_file\n" >> GeCo3_Output/Archaea/PT/ARCHAEA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Archaea/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Archaea/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/PT/geco_ARCHAEA.txt
                rm GeCo3_Output/Archaea/PT/TO_COMPRESS*
            done
        fi
    done
}
#
# ################## FUNGI ##################
#
function FUNGI_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Fungi Genomes"

    if [[ ! -d "GeCo3_Output/Fungi" ]]; then
        mkdir GeCo3_Output/Fungi
    fi

    if [[ ! -d "GeCo3_Output/Fungi/PT" ]]; then
        mkdir GeCo3_Output/Fungi/PT
    fi

    if [[ -f "GeCo3_Output/Fungi/PT/FUNGI_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Fungi/PT/FUNGI_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Fungi/PT-Fungi/*
    #for file in References/NCBI-Fungi/PT-Fungi/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "fungi_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Fungi/PT/FUNGI_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Fungi/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Fungi/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/PT/geco_FUNGI.txt
                rm GeCo3_Output/Fungi/PT/TO_COMPRESS*
            done
        fi
    done
}
#
# ################## PLANT ##################
#
function PLANT_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plant Genomes"

    if [[ ! -d "GeCo3_Output/Plant" ]]; then
        mkdir GeCo3_Output/Plant
    fi

    if [[ ! -d "GeCo3_Output/Plant/PT" ]]; then
        mkdir GeCo3_Output/Plant/PT
    fi

    if [[ -f "GeCo3_Output/Plant/PT/PLANT_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Plant/PT/PLANT_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Plant/PT-Plant/*
    #for file in References/NCBI-Plant/PT-Plant/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "plant_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Plant/PT/PLANT_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plant/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Plant/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/PT/geco_PLANT.txt
                rm GeCo3_Output/Plant/PT/TO_COMPRESS*
            done
        fi
    done
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Protozoa Genomes"

    if [[ ! -d "GeCo3_Output/Protozoa" ]]; then
        mkdir GeCo3_Output/Protozoa
    fi

    if [[ ! -d "GeCo3_Output/Protozoa/PT" ]]; then
        mkdir GeCo3_Output/Protozoa/PT
    fi

    if [[ -f "GeCo3_Output/Protozoa/PT/PROTOZOA_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Protozoa/PT/PROTOZOA_ID.txt
    fi

    #for file in /media/alexloure/T7Touch/NCBI-Protozoa/PT-Protozoa/*
    for file in References/NCBI-Protozoa/PT-Protozoa/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Protozoa/PT/PROTOZOA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Protozoa/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Protozoa/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/PT/geco_PROTOZOA.txt
                rm GeCo3_Output/Protozoa/PT/TO_COMPRESS*
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
    mkdir $path/PT-plastid
    awk '/>/ {filename="/media/alexloure/T7Touch/NCBI-Plastid/PT-plastid/"NR".fna"} {print > filename}' $path/DB-plastid.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plastid Genomes"

    if [[ ! -d "GeCo3_Output/Plastid" ]]; then
        mkdir GeCo3_Output/Plastid
    fi

    if [[ ! -d "GeCo3_Output/Plastid/PT" ]]; then
        mkdir GeCo3_Output/Plastid/PT
    fi

    if [[ -f "GeCo3_Output/Plastid/PT/PLASTID_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Plastid/PT/PLASTID_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Plastid/PT-Plastid/*
    #for file in References/NCBI-Plastid/PT-Plastid/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "plastid_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Plastid/PT/PLASTID_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plastid/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Plastid/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plastid/PT/geco_PLASTID.txt
                rm GeCo3_Output/Plastid/PT/TO_COMPRESS*
            done
        fi
    done
    rm References/NCBI-Plastid/DB-plastid.fa
}
#
# ################## MITOCHONDRIAL ##################
#
function MITOCHONDRIAL_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start parsing Mitochondrial Genomes"
    ulimit -n 1000000
    path="/media/alexloure/T7Touch/NCBI-Mitochondrial"
    gunzip -k $path/DB-mitochondrion.fa.gz
    mkdir $path/PT-mitochondrion
    awk '/>/ {filename="/media/alexloure/T7Touch/NCBI-Mitochondrial/PT-mitochondrion/"NR".fna"} {print > filename}' $path/DB-mitochondrion.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

    if [[ ! -d "GeCo3_Output/Mitochondrial" ]]; then
        mkdir GeCo3_Output/Mitochondrial
    fi

    if [[ ! -d "GeCo3_Output/Mitochondrial/PT" ]]; then
        mkdir GeCo3_Output/Mitochondrial/PT
    fi

    if [[ -f "GeCo3_Output/Mitochondrial/PT/MITOCHONDRIAL_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Mitochondrial/PT/MITOCHONDRIAL_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Mitochondrial/PT-Mitochondrial/*
    #for file in References/NCBI-Mitochondrial/PT-Mitochondrial/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Mitochondrial/PT/MITOCHONDRIAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Mitochondrial/PT/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Mitochondrial/PT/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Mitochondrial/PT/geco_MITOCHONDRIAL.txt
                rm GeCo3_Output/Mitochondrial/PT/TO_COMPRESS*
            done
        fi
    done
    rm References/NCBI-Mitochondrial/DB-mitochondrion.fa
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Viral DB"
  ./src/compression_analysis_csv.sh Virus VIRAL geco3_Viral.csv PT
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_GENOME_COMPRESSION;
    echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Bacteria DB"
  ./src/compression_analysis_csv.sh Bacteria BACTERIA geco3_Bacteria.csv PT
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Archeae DB"
  ./src/compression_analysis_csv.sh Archaea ARCHAEA geco3_Archaea.csv PT
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Protozoa DB"
  ./src/compression_analysis_csv.sh Protozoa PROTOZOA geco3_Protozoa.csv PT
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Fungi DB"
  ./src/compression_analysis_csv.sh Fungi FUNGI geco3_Fungi.csv PT
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plant DB"
  ./src/compression_analysis_csv.sh Plant PLANT geco3_Plant.csv PT
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Mitochondrial DB"
  ./src/compression_analysis_csv.sh Mitochondrial MITOCHONDRIAL geco3_Mitochondrial.csv PT
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plastid DB"
  ./src/compression_analysis_csv.sh Plastid PLASTID geco3_Plastid.csv PT
  fi
#
# ==============================================================================
#