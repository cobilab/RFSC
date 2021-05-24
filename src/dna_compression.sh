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

    if [[ ! -d "GeCo3_Output/Virus" ]]; then
        mkdir GeCo3_Output/Virus
    fi

    if [[ ! -d "GeCo3_Output/Virus/NM" ]]; then
        mkdir GeCo3_Output/Virus/NM
    fi

    if [[ -f "GeCo3_Output/Virus/NM/VIRAL_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Virus/NM/VIRAL_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Virus/NM-viral/*
    for file in References/NCBI-Virus/NM-viral/*
    do
        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "viral_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";
        
        else
            printf "$out_file\n" >> GeCo3_Output/Virus/NM/VIRAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Virus/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Virus/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/NM/geco_VIRAL.txt
                rm GeCo3_Output/Virus/NM/TO_COMPRESS*
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

    if [[ ! -d "GeCo3_Output/Bacteria/NM" ]]; then
        mkdir GeCo3_Output/Bacteria/NM
    fi

    if [[ -f "GeCo3_Output/Bacteria/NM/BACTERIA_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Bacteria/NM/BACTERIA_ID.txt
    fi

    for file in /media/alexloure/T7Touch/NCBI-Bacteria/NM-bacteria/*
    # for file in References/NCBI-Bacteria/NM-bacteria/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Bacteria/NM/BACTERIA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Bacteria/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Bacteria/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/NM/geco_BACTERIA.txt
                rm GeCo3_Output/Bacteria/NM/TO_COMPRESS*
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

    if [[ ! -d "GeCo3_Output/Archaea/NM" ]]; then
        mkdir GeCo3_Output/Archaea/NM
    fi

    if [[ -f "GeCo3_Output/Archaea/NM/ARCHAEA_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Archaea/NM/ARCHAEA_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Archaea/NM-archaea/*
    for file in References/NCBI-Archaea/NM-archaea/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "archaea_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        
        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";
        
        else
            printf "$out_file\n" >> GeCo3_Output/Archaea/NM/ARCHAEA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Archaea/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Archaea/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/NM/geco_ARCHAEA.txt
                rm GeCo3_Output/Archaea/NM/TO_COMPRESS*
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

    if [[ ! -d "GeCo3_Output/Fungi/NM" ]]; then
        mkdir GeCo3_Output/Fungi/NM
    fi

    if [[ -f "GeCo3_Output/Fungi/NM/FUNGI_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Fungi/NM/FUNGI_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Fungi/NM-fungi/*
    for file in References/NCBI-Fungi/NM-fungi/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "fungi_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Fungi/NM/FUNGI_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Fungi/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Fungi/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/NM/geco_FUNGI.txt
                rm GeCo3_Output/Fungi/NM/TO_COMPRESS*
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

    if [[ ! -d "GeCo3_Output/Plant/NM" ]]; then
        mkdir GeCo3_Output/Plant/NM
    fi

    if [[ -f "GeCo3_Output/Plant/NM/PLANT_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Plant/NM/PLANT_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Plant/NM-plant/*
    for file in References/NCBI-Plant/NM-plant/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "plant_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Plant/NM/PLANT_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plant/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Plant/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/NM/geco_PLANT.txt
                rm GeCo3_Output/Plant/NM/TO_COMPRESS*
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

    if [[ ! -d "GeCo3_Output/Protozoa/NM" ]]; then
        mkdir GeCo3_Output/Protozoa/NM
    fi

    if [[ -f "GeCo3_Output/Protozoa/NM/PROTOZOA_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Protozoa/NM/PROTOZOA_ID.txt
    fi

    #for file in /media/alexloure/T7Touch/NCBI-Protozoa/NM-protozoa/*
    for file in References/NCBI-Protozoa/NM-protozoa/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Protozoa/NM/PROTOZOA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Protozoa/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Protozoa/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/NM/geco_PROTOZOA.txt
                rm GeCo3_Output/Protozoa/NM/TO_COMPRESS*
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

    if [[ ! -d "GeCo3_Output/Plastid" ]]; then
        mkdir GeCo3_Output/Plastid
    fi

    if [[ ! -d "GeCo3_Output/Plastid/NM" ]]; then
        mkdir GeCo3_Output/Plastid/NM
    fi

    if [[ -f "GeCo3_Output/Plastid/NM/PLASTID_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Plastid/NM/PLASTID_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Plastid/NM-plastid/*
    for file in References/NCBI-Plastid/NM-plastid/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "plastid_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Plastid/NM/PLASTID_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plastid/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Plastid/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plastid/NM/geco_PLASTID.txt
                rm GeCo3_Output/Plastid/NM/TO_COMPRESS*
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
    mkdir $path/NM-mitochondrion
    awk '/>/ {filename="/media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/"NR".fna"} {print > filename}' $path/DB-mitochondrion.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

    if [[ ! -d "GeCo3_Output/Mitochondrial" ]]; then
        mkdir GeCo3_Output/Mitochondrial
    fi

    if [[ ! -d "GeCo3_Output/Mitochondrial/NM" ]]; then
        mkdir GeCo3_Output/Mitochondrial/NM
    fi

    if [[ -f "GeCo3_Output/Mitochondrial/NM/MITOCHONDRIAL_ID.txt" ]]; then
        readarray -t already_processed < GeCo3_Output/Mitochondrial/NM/MITOCHONDRIAL_ID.txt
    fi

    # for file in /media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/*
    for file in References/NCBI-Mitochondrial/NM-mitochondrion/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            printf "$out_file\n" >> GeCo3_Output/Mitochondrial/NM/MITOCHONDRIAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Mitochondrial/NM/TO_COMPRESS
                GeCo3 -v -l $l GeCo3_Output/Mitochondrial/NM/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Mitochondrial/NM/geco_MITOCHONDRIAL.txt
                rm GeCo3_Output/Mitochondrial/NM/TO_COMPRESS*
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
  ./src/compression_analysis_csv.sh Virus VIRAL geco3_Viral.csv NM
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_GENOME_COMPRESSION;
    echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Bacteria DB"
  ./src/compression_analysis_csv.sh Bacteria BACTERIA geco3_Bacteria.csv NM
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Archeae DB"
  ./src/compression_analysis_csv.sh Archaea ARCHAEA geco3_Archaea.csv NM
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Protozoa DB"
  ./src/compression_analysis_csv.sh Protozoa PROTOZOA geco3_Protozoa.csv NM
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Fungi DB"
  ./src/compression_analysis_csv.sh Fungi FUNGI geco3_Fungi.csv NM
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plant DB"
  ./src/compression_analysis_csv.sh Plant PLANT geco3_Plant.csv NM
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Mitochondrial DB"
  ./src/compression_analysis_csv.sh Mitochondrial MITOCHONDRIAL geco3_Mitochondrial.csv NM
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plastid DB"
  ./src/compression_analysis_csv.sh Plastid PLASTID geco3_Plastid.csv NM
  fi
#
# ==============================================================================
#