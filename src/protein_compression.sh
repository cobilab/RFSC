#!/bin/bash
#
testing_levels=(1 2 3 4 7);
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

    if [[ ! -d "Analysis/AC/Virus" ]]; then
        mkdir Analysis/AC/Virus
    fi

    if [[ ! -d "Analysis/AC/Virus" ]]; then
        mkdir Analysis/AC/Virus
    fi

    if [[ -f "Analysis/AC/Virus/VIRAL_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Virus/VIRAL_ID.txt
    fi

    for file in References/NCBI-Virus/PT-Virus/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "viral_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Virus/PT-Virus/*.co
            printf "$out_file\n" >> Analysis/AC/Virus/VIRAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Virus/ac_VIRAL.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Virus/check_virus.txt
            done
        fi
    done
}
#
################## BACTERIA ##################
#
function BACTERIA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Bacteria Genomes"

    if [[ ! -d "Analysis/AC/Bacteria" ]]; then
        mkdir Analysis/AC/Bacteria
    fi

    if [[ ! -d "Analysis/AC/Bacteria" ]]; then
        mkdir Analysis/AC/Bacteria
    fi

    if [[ -f "Analysis/AC/Bacteria/BACTERIA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Bacteria/BACTERIA_ID.txt
    fi

    for file in References/NCBI-Bacteria/PT-Bacteria/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Bacteria/PT-Bacteria/*.co
            printf "$out_file\n" >> Analysis/AC/Bacteria/BACTERIA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Bacteria/ac_BACTERIA.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Bacteria/check_bacteria.txt
            done
        fi
    done   
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Archaea Genomes"

    if [[ ! -d "Analysis/AC/Archaea" ]]; then
        mkdir Analysis/AC/Archaea
    fi

    if [[ ! -d "Analysis/AC/Archaea" ]]; then
        mkdir Analysis/AC/Archaea
    fi

    if [[ -f "Analysis/AC/Archaea/ARCHAEA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Archaea/ARCHAEA_ID.txt
    fi

    for file in References/NCBI-Archaea/PT-Archaea/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "archaea_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        
        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";
        
        else
            rm References/NCBI-Archaea/PT-Archaea/*.co
            printf "$out_file\n" >> Analysis/AC/Archaea/ARCHAEA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Archaea/ac_ARCHAEA.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Archaea/check_archaea.txt
            done
        fi
    done
}
#
# ################## FUNGI ##################
#
function FUNGI_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Fungi Genomes"

    if [[ ! -d "Analysis/AC/Fungi" ]]; then
        mkdir Analysis/AC/Fungi
    fi

    if [[ ! -d "Analysis/AC/Fungi" ]]; then
        mkdir Analysis/AC/Fungi
    fi

    if [[ -f "Analysis/AC/Fungi/FUNGI_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Fungi/FUNGI_ID.txt
    fi

    for file in References/NCBI-Fungi/PT-Fungi/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "fungi_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Fungi/PT-Fungi/*.co
            printf "$out_file\n" >> Analysis/AC/Fungi/FUNGI_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Fungi/ac_FUNGI.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Fungi/check_fungi.txt
            done
        fi
    done
}
#
# ################## PLANT ##################
#
function PLANT_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plant Genomes"

    if [[ ! -d "Analysis/AC/Plant" ]]; then
        mkdir Analysis/AC/Plant
    fi

    if [[ ! -d "Analysis/AC/Plant" ]]; then
        mkdir Analysis/AC/Plant
    fi

    if [[ -f "Analysis/AC/Plant/PLANT_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Plant/PLANT_ID.txt
    fi

    for file in References/NCBI-Plant/PT-Plant/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "plant_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Plant/PT-Plant/*.co
            printf "$out_file\n" >> Analysis/AC/Plant/PLANT_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Plant/ac_PLANT.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Plant/check_plant.txt
            done
        fi
    done
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Protozoa Genomes"

    if [[ ! -d "Analysis/AC/Protozoa" ]]; then
        mkdir Analysis/AC/Protozoa
    fi

    if [[ ! -d "Analysis/AC/Protozoa" ]]; then
        mkdir Analysis/AC/Protozoa
    fi

    if [[ -f "Analysis/AC/Protozoa/PROTOZOA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Protozoa/PROTOZOA_ID.txt
    fi

    for file in References/NCBI-Protozoa/PT-Protozoa/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Protozoa/PT-Protozoa/*.co
            printf "$out_file\n" >> Analysis/AC/Protozoa/PROTOZOA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Protozoa/ac_PROTOZOA.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Protozoa/check_protozoa.txt
            done
        fi
    done
}
#
################## PLASTID ##################
#
function PLASTID_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plastid Genomes"

    if [[ ! -d "Analysis/AC/Plastid" ]]; then
        mkdir Analysis/AC/Plastid
    fi

    if [[ ! -d "Analysis/AC/Plastid" ]]; then
        mkdir Analysis/AC/Plastid
    fi

    if [[ -f "Analysis/AC/Plastid/PLASTID_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Plastid/PLASTID_ID.txt
    fi

    for file in References/NCBI-Plastid/PT-Plastid/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "plastid_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Plastid/PT-Plastid/*.co
            printf "$out_file\n" >> Analysis/AC/Plastid/PLASTID_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Plastid/ac_PLASTID.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Plastid/check_plastid.txt
            done
        fi
    done
}
#
# ################## MITOCHONDRIAL ##################
#
function MITOCHONDRIAL_GENOME_COMPRESSION () {

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

    if [[ ! -d "Analysis/AC/Mitochondrial" ]]; then
        mkdir Analysis/AC/Mitochondrial
    fi

    if [[ ! -d "Analysis/AC/Mitochondrial" ]]; then
        mkdir Analysis/AC/Mitochondrial
    fi

    if [[ -f "Analysis/AC/Mitochondrial/MITOCHONDRIAL_ID.txt" ]]; then
        readarray -t already_processed < Analysis/AC/Mitochondrial/MITOCHONDRIAL_ID.txt
    fi

    for file in References/NCBI-Mitochondrial/PT-Mitochondrial/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Mitochondrial/PT-Mitochondrial/*.co
            printf "$out_file\n" >> Analysis/AC/Mitochondrial/MITOCHONDRIAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                AC -v -l $l $file | sed '3,4d' | sed 'N;s/\n/ /' >> Analysis/AC/Mitochondrial/ac_MITOCHONDRIAL.txt
                rm $file.co

                echo "$out_file $l\n" >> Analysis/AC/Mitochondrial/check_mito.txt
            done
        fi
    done
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Viral DB"
  ./src/compression_analysis_csv.sh Virus VIRAL ac2_Viral.csv PT
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_GENOME_COMPRESSION;
    echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Bacteria DB"
  ./src/compression_analysis_csv.sh Bacteria BACTERIA ac2_Bacteria.csv PT
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Archeae DB"
  ./src/compression_analysis_csv.sh Archaea ARCHAEA ac2_Archaea.csv PT
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Protozoa DB"
  ./src/compression_analysis_csv.sh Protozoa PROTOZOA ac2_Protozoa.csv PT
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Fungi DB"
  ./src/compression_analysis_csv.sh Fungi FUNGI ac2_Fungi.csv PT
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plant DB"
  ./src/compression_analysis_csv.sh Plant PLANT ac2_Plant.csv PT
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Mitochondrial DB"
  ./src/compression_analysis_csv.sh Mitochondrial MITOCHONDRIAL ac2_Mitochondrial.csv PT
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_GENOME_COMPRESSION;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plastid DB"
  ./src/compression_analysis_csv.sh Plastid PLASTID ac2_Plastid.csv PT
  fi
#
# ==============================================================================
#