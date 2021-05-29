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
  echo "Usage: ./GCcontent_DB_analysis.sh [options]                                     "
  echo "                                                                                "
  echo "It analyses the percentage of GC-Content in each sequece                        "
  echo "of the choosen NCBI database:                                                   "
  echo "                                                                                "
  echo "   -h,   --help                   Show this help message,                       "
  echo "                                                                                "
  echo "   -vi,  --viral                  GC-Content viral DB,                          "
  echo "   -ba,  --bacteria               GC-Content bacteria DB,                       "
  echo "   -ar,  --archaea                GC-Content archaea DB,                        "
  echo "   -pr,  --protozoa               GC-Content protozoa DB,                       "
  echo "   -fu,  --fungi                  GC-Content fungi DB,                          "
  echo "   -pl,  --plant                  GC-Content plant DB,                          "
  echo "   -mi,  --mitochondrion          GC-Content mitochondrion DB,                  "
  echo "   -ps,  --plastid                GC-Content plastid DB.                        "
  echo "                                                                                "
  echo "Example: ./GCcontent_DB_analysis.sh --viral --fungi --plant                     "
  echo "                                                                                "
  exit 1
  fi
#
# ==============================================================================
#
################## VIRUS ##################
#
function VIRUS_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Viral GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Virus" ]]; then
        mkdir Analysis/GCcontent/Virus
    fi

    if [[ -f "Analysis/GCcontent/Virus/VIRAL_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Virus/VIRAL_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Virus/VIRAL_ID.txt
            ./src/GC_analysis.sh $file ZIPPED >> Analysis/GCcontent/Virus/GC_content_VIRAL.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Viral GC-Content analysis!\033[0m";
}
#
################## BACTERIA ##################
#
function BACTERIA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Bacteria GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Bacteria" ]]; then
        mkdir Analysis/GCcontent/Bacteria
    fi

    if [[ -f "Analysis/GCcontent/Bacteria/NM/BACTERIA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Bacteria/NM/BACTERIA_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Bacteria/BACTERIA_ID.txt
            ./src/GC_analysis.sh $file ZIPPED >> Analysis/GCcontent/Bacteria/GC_content_BACTERIA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Bacteria GC-Content analysis!\033[0m";
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Archaea GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Archaea" ]]; then
        mkdir Analysis/GCcontent/Archaea
    fi

    if [[ -f "Analysis/GCcontent/Archaea/ARCHAEA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Archaea/ARCHAEA_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Archaea/ARCHAEA_ID.txt
            ./src/GC_analysis.sh $file ZIPPED >> Analysis/GCcontent/Archaea/GC_content_ARCHAEA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Archaea GC-Content analysis!\033[0m";
}
#
# ################## FUNGI ##################
#
function FUNGI_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Fungi GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Fungi" ]]; then
        mkdir Analysis/GCcontent/Fungi
    fi

    if [[ -f "Analysis/GCcontent/Fungi/FUNGI_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Fungi/FUNGI_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Fungi/FUNGI_ID.txt
            ./src/GC_analysis.sh $file ZIPPED >> Analysis/GCcontent/Fungi/GC_content_FUNGI.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Fungi GC-Content analysis!\033[0m";
}
#
# ################## PLANT ##################
#
function PLANT_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Plant GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Plant" ]]; then
        mkdir Analysis/GCcontent/Plant
    fi

    if [[ -f "Analysis/GCcontent/Plant/PLANT_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Plant/PLANT_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Plant/PLANT_ID.txt
            ./src/GC_analysis.sh $file ZIPPED >> Analysis/GCcontent/Plant/GC_content_PLANT.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Plant GC-Content analysis!\033[0m";
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Protozoa GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Protozoa" ]]; then
        mkdir Analysis/GCcontent/Protozoa
    fi

    if [[ -f "Analysis/GCcontent/Protozoa/PROTOZOA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Protozoa/PROTOZOA_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Protozoa/PROTOZOA_ID.txt
            ./src/GC_analysis.sh $file ZIPPED >> Analysis/GCcontent/Protozoa/GC_content_PROTOZOA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Protozoa GC-Content analysis!\033[0m";
}
#
################## PLASTID ##################
#
function PLASTID_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Plastid GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Plastid" ]]; then
        mkdir Analysis/GCcontent/Plastid
    fi

    if [[ -f "Analysis/GCcontent/Plastid/PLASTID_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Plastid/PLASTID_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Plastid/PLASTID_ID.txt
            ./src/GC_analysis.sh $file UNZIPPED >> Analysis/GCcontent/Plastid/GC_content_PLASTID.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Plastid GC-Content analysis!\033[0m";
}
#
# ################## MITOCHONDRIAL ##################
#
function MITOCHONDRIAL_GENOME_COMPRESSION () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Mitochondrial GC-Content Analysis"

    if [[ ! -d "Analysis/GCcontent" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/GCcontent was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/GCcontent/Mitochondrial" ]]; then
        mkdir Analysis/GCcontent/Mitochondrial
    fi

    if [[ -f "Analysis/GCcontent/Mitochondrial/MITOCHONDRIAL_ID.txt" ]]; then
        readarray -t already_processed < Analysis/GCcontent/Mitochondrial/MITOCHONDRIAL_ID.txt
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
            printf "$out_file\n" >> Analysis/GCcontent/Mitochondrial/MITOCHONDRIAL_ID.txt
            ./src/GC_analysis.sh $file UNZIPPED >> Analysis/GCcontent/Mitochondrial/GC_content_MITOCHONDRIAL.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Mitochondrial GC-Content analysis!\033[0m";
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_GENOME_COMPRESSION;
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_GENOME_COMPRESSION;
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_GENOME_COMPRESSION;
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_GENOME_COMPRESSION;
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_GENOME_COMPRESSION;
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_GENOME_COMPRESSION;
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_GENOME_COMPRESSION;
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_GENOME_COMPRESSION;
  fi
#
# ==============================================================================
#