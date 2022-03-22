#!/bin/bash
#
#testing_levels=(1 2 3 4 7 12 14);
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

    if [[ ! -d "Analysis/GeCo/Virus" ]]; then
        mkdir Analysis/GeCo/Virus
    fi

    if [[ ! -d "Analysis/GeCo/Virus" ]]; then
        mkdir Analysis/GeCo/Virus
    fi

    if [[ -f "Analysis/GeCo/Virus/VIRAL_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Virus/VIRAL_ID.txt
    fi

    for file in References/NCBI-Virus/NM-viral/*
    do
        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "viral_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";
        
        else
            rm References/NCBI-Virus/NM-viral/*.co
            printf "$out_file\n" >> Analysis/GeCo/Virus/VIRAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Virus/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Virus/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Virus/geco_VIRAL.txt
                rm Analysis/GeCo/Virus/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Virus/check_virus.txt
            done
        fi
    done
}
#
################## BACTERIA ##################
#
function BACTERIA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Bacteria Genomes"

    if [[ ! -d "Analysis/GeCo/Bacteria" ]]; then
        mkdir Analysis/GeCo/Bacteria
    fi

    if [[ ! -d "Analysis/GeCo/Bacteria" ]]; then
        mkdir Analysis/GeCo/Bacteria
    fi

    if [[ -f "Analysis/GeCo/Bacteria/BACTERIA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Bacteria/BACTERIA_ID.txt
    fi

    for file in References/NCBI-Bacteria/NM-bacteria/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Bacteria/NM-bacteria/*.co
            printf "$out_file\n" >> Analysis/GeCo/Bacteria/BACTERIA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Bacteria/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Bacteria/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Bacteria/geco_BACTERIA.txt
                rm Analysis/GeCo/Bacteria/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Bacteria/check_bacteria.txt
            done
        fi
    done   
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Archaea Genomes"

    if [[ ! -d "Analysis/GeCo/Archaea" ]]; then
        mkdir Analysis/GeCo/Archaea
    fi

    if [[ ! -d "Analysis/GeCo/Archaea" ]]; then
        mkdir Analysis/GeCo/Archaea
    fi

    if [[ -f "Analysis/GeCo/Archaea/ARCHAEA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Archaea/ARCHAEA_ID.txt
    fi

    for file in References/NCBI-Archaea/NM-archaea/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "archaea_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
        
        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";
        
        else
            rm References/NCBI-Archaea/NM-archaea/*.co
            printf "$out_file\n" >> Analysis/GeCo/Archaea/ARCHAEA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Archaea/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Archaea/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Archaea/geco_ARCHAEA.txt
                rm Analysis/GeCo/Archaea/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Archaea/check_archaea.txt
            done
        fi
    done
}
#
# ################## FUNGI ##################
#
function FUNGI_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Fungi Genomes"

    if [[ ! -d "Analysis/GeCo/Fungi" ]]; then
        mkdir Analysis/GeCo/Fungi
    fi

    if [[ ! -d "Analysis/GeCo/Fungi" ]]; then
        mkdir Analysis/GeCo/Fungi
    fi

    if [[ -f "Analysis/GeCo/Fungi/FUNGI_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Fungi/FUNGI_ID.txt
    fi

    for file in References/NCBI-Fungi/NM-fungi/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "fungi_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Fungi/NM-fungi/*.co
            printf "$out_file\n" >> Analysis/GeCo/Fungi/FUNGI_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Fungi/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Fungi/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Fungi/geco_FUNGI.txt
                rm Analysis/GeCo/Fungi/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Fungi/check_fungi.txt
            done
        fi
    done
}
#
# ################## PLANT ##################
#
function PLANT_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plant Genomes"

    if [[ ! -d "Analysis/GeCo/Plant" ]]; then
        mkdir Analysis/GeCo/Plant
    fi

    if [[ ! -d "Analysis/GeCo/Plant" ]]; then
        mkdir Analysis/GeCo/Plant
    fi

    if [[ -f "Analysis/GeCo/Plant/PLANT_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Plant/PLANT_ID.txt
    fi

    for file in References/NCBI-Plant/NM-plant/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "plant_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Plant/NM-plant/*.co
            printf "$out_file\n" >> Analysis/GeCo/Plant/PLANT_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Plant/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Plant/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Plant/geco_PLANT.txt
                rm Analysis/GeCo/Plant/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Plant/check_plant.txt
            done
        fi
    done
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Protozoa Genomes"

    if [[ ! -d "Analysis/GeCo/Protozoa" ]]; then
        mkdir Analysis/GeCo/Protozoa
    fi

    if [[ ! -d "Analysis/GeCo/Protozoa" ]]; then
        mkdir Analysis/GeCo/Protozoa
    fi

    if [[ -f "Analysis/GeCo/Protozoa/PROTOZOA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Protozoa/PROTOZOA_ID.txt
    fi

    for file in References/NCBI-Protozoa/NM-protozoa/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna.gz"}"

        if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Protozoa/NM-protozoa/*.co
            printf "$out_file\n" >> Analysis/GeCo/Protozoa/PROTOZOA_ID.txt

            for l in ${testing_levels[@]}; 
            do
                zcat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Protozoa/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Protozoa/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Protozoa/geco_PROTOZOA.txt
                rm Analysis/GeCo/Protozoa/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Protozoa/check_protozoa.txt
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
    path="References/NCBI-Plastid"
    gunzip -k $path/DB-plastid.fa.gz
    mkdir $path-plastid
    awk '/>/ {filename="References/NCBI-Plastid-plastid/"NR".fna"} {print > filename}' $path/DB-plastid.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plastid Genomes"

    if [[ ! -d "Analysis/GeCo/Plastid" ]]; then
        mkdir Analysis/GeCo/Plastid
    fi

    if [[ ! -d "Analysis/GeCo/Plastid" ]]; then
        mkdir Analysis/GeCo/Plastid
    fi

    if [[ -f "Analysis/GeCo/Plastid/PLASTID_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Plastid/PLASTID_ID.txt
    fi

    for file in References/NCBI-Plastid/NM-plastid/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "plastid_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Plastid/NM-plastid/*.co
            printf "$out_file\n" >> Analysis/GeCo/Plastid/PLASTID_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Plastid/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Plastid/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Plastid/geco_PLASTID.txt
                rm Analysis/GeCo/Plastid/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Plastid/check_plastid.txt
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
    path="References/NCBI-Mitochondrial"
    gunzip -k $path/DB-mitochondrion.fa.gz
    mkdir $path-mitochondrion
    awk '/>/ {filename="References/NCBI-Mitochondrial-mitochondrion/"NR".fna"} {print > filename}' $path/DB-mitochondrion.fa

    echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

    if [[ ! -d "Analysis/GeCo/Mitochondrial" ]]; then
        mkdir Analysis/GeCo/Mitochondrial
    fi

    if [[ ! -d "Analysis/GeCo/Mitochondrial" ]]; then
        mkdir Analysis/GeCo/Mitochondrial
    fi

    if [[ -f "Analysis/GeCo/Mitochondrial/MITOCHONDRIAL_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GeCo/Mitochondrial/MITOCHONDRIAL_ID.txt
    fi

    for file in References/NCBI-Mitochondrial/NM-mitochondrion/*
    do

        out_file=$(basename $file);
        out_file="${out_file%".fna"}"

        if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";

        elif [[ " ${already_processed[@]} " =~ " ${out_file} " ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m File ($out_file) already was processed";

        else
            rm References/NCBI-Mitochondrial/NM-mitochondrion/*.co
            printf "$out_file\n" >> Analysis/GeCo/Mitochondrial/MITOCHONDRIAL_ID.txt

            for l in ${testing_levels[@]}; 
            do
                cat $file | grep -v ">" | tr -d -c "ACGT" > Analysis/GeCo/Mitochondrial/TO_COMPRESS
                GeCo3 -v -l $l Analysis/GeCo/Mitochondrial/TO_COMPRESS | sed '1,6d' | sed '2d' >> Analysis/GeCo/Mitochondrial/geco_MITOCHONDRIAL.txt
                rm Analysis/GeCo/Mitochondrial/TO_COMPRESS*

                echo "$out_file $l\n" >> Analysis/GeCo/Mitochondrial/check_mito.txt
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