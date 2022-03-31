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
  echo "Usage: ./length_DB_analysis.sh [options]                                        "
  echo "                                                                                "
  echo "It analyses the Length sequences the choosen NCBI database:                     "
  echo "                                                                                "
  echo "   -h,   --help                   Show this help message,                       "
  echo "                                                                                "
  echo "   -vi,  --viral                  Length sequences in viral DB,                 "
  echo "   -ba,  --bacteria               Length sequences in bacteria DB,              "
  echo "   -ar,  --archaea                Length sequences in archaea DB,               "
  echo "   -pr,  --protozoa               Length sequences in protozoa DB,              "
  echo "   -fu,  --fungi                  Length sequences in fungi DB,                 "
  echo "   -pl,  --plant                  Length sequences in plant DB,                 "
  echo "   -mi,  --mitochondrion          Length sequences in mitochondrion DB,         "
  echo "   -ps,  --plastid                Length sequences in plastid DB.               "
  echo "                                                                                "
  echo "Example: ./length_DB_analysis.sh --viral --fungi --plant                        "
  echo "                                                                                "
  exit 1
  fi
#
# ==============================================================================
#
################## VIRUS ##################
#
function VIRUS_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Viral DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Virus" ]]; then
        mkdir Analysis/LengthSeq/Virus
    fi

    if [[ -f "Analysis/LengthSeq/Virus/VIRAL_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Virus/VIRAL_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Virus/VIRAL_DNA_ID.txt
            ./src/length_sequence.sh $file DNA ZIPPED >> Analysis/LengthSeq/Virus/Length_DNA_VIRAL.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Viral DNA Length analysis!\033[0m";
}
function VIRUS_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Viral AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Virus" ]]; then
        mkdir Analysis/LengthSeq/Virus
    fi

    if [[ -f "Analysis/LengthSeq/Virus/VIRAL_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Virus/VIRAL_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Virus/VIRAL_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Virus/Length_AA_VIRAL.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Viral AA Length analysis!\033[0m";
}
#
################## BACTERIA ##################
#
function BACTERIA_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Bacteria DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Bacteria" ]]; then
        mkdir Analysis/LengthSeq/Bacteria
    fi

    if [[ -f "Analysis/LengthSeq/Bacteria/NM/BACTERIA_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Bacteria/NM/BACTERIA_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Bacteria/BACTERIA_DNA_ID.txt
            ./src/length_sequence.sh $file DNA ZIPPED >> Analysis/LengthSeq/Bacteria/Length_DNA_BACTERIA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Bacteria DNA Length analysis!\033[0m";
}
function BACTERIA_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Bacteria AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Bacteria" ]]; then
        mkdir Analysis/LengthSeq/Bacteria
    fi

    if [[ -f "Analysis/LengthSeq/Bacteria/NM/BACTERIA_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Bacteria/NM/BACTERIA_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Bacteria/BACTERIA_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Bacteria/Length_AA_BACTERIA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Bacteria AA Length analysis!\033[0m";
}
#
# ################## ARCHAEA ##################
#
function ARCHAEA_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Archaea DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Archaea" ]]; then
        mkdir Analysis/LengthSeq/Archaea
    fi

    if [[ -f "Analysis/LengthSeq/Archaea/ARCHAEA_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Archaea/ARCHAEA_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Archaea/ARCHAEA_DNA_ID.txt
            ./src/length_sequence.sh $file DNA ZIPPED >> Analysis/LengthSeq/Archaea/Length_DNA_ARCHAEA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Archaea DNA Length analysis!\033[0m";
}
function ARCHAEA_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Archaea AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Archaea" ]]; then
        mkdir Analysis/LengthSeq/Archaea
    fi

    if [[ -f "Analysis/LengthSeq/Archaea/ARCHAEA_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Archaea/ARCHAEA_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Archaea/ARCHAEA_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Archaea/Length_AA_ARCHAEA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Archaea AA Length analysis!\033[0m";
}
#
# ################## FUNGI ##################
#
function FUNGI_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Fungi DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Fungi" ]]; then
        mkdir Analysis/LengthSeq/Fungi
    fi

    if [[ -f "Analysis/LengthSeq/Fungi/FUNGI_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Fungi/FUNGI_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Fungi/FUNGI_DNA_ID.txt
            ./src/length_sequence.sh $file DNA ZIPPED >> Analysis/LengthSeq/Fungi/Length_DNA_FUNGI.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Fungi DNA Length analysis!\033[0m";
}
function FUNGI_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Fungi AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Fungi" ]]; then
        mkdir Analysis/LengthSeq/Fungi
    fi

    if [[ -f "Analysis/LengthSeq/Fungi/FUNGI_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Fungi/FUNGI_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Fungi/FUNGI_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Fungi/Length_AA_FUNGI.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Fungi AA Length analysis!\033[0m";
}
#
# ################## PLANT ##################
#
function PLANT_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Plant DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Plant" ]]; then
        mkdir Analysis/LengthSeq/Plant
    fi

    if [[ -f "Analysis/LengthSeq/Plant/PLANT_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Plant/PLANT_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Plant/PLANT_DNA_ID.txt
            ./src/length_sequence.sh $file DNA ZIPPED >> Analysis/LengthSeq/Plant/Length_DNA_PLANT.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Plant DNA Length analysis!\033[0m";
}
function PLANT_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Plant AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Plant" ]]; then
        mkdir Analysis/LengthSeq/Plant
    fi

    if [[ -f "Analysis/LengthSeq/Plant/PLANT_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Plant/PLANT_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Plant/PLANT_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Plant/Length_AA_PLANT.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Plant AA Length analysis!\033[0m";
}
#
# ################## PROTOZOA ##################
#
function PROTOZOA_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Protozoa DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Protozoa" ]]; then
        mkdir Analysis/LengthSeq/Protozoa
    fi

    if [[ -f "Analysis/LengthSeq/Protozoa/PROTOZOA_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Protozoa/PROTOZOA_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Protozoa/PROTOZOA_DNA_ID.txt
            ./src/length_sequence.sh $file DNA ZIPPED >> Analysis/LengthSeq/Protozoa/Length_DNA_PROTOZOA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Protozoa DNA Length analysis!\033[0m";
}
function PROTOZOA_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Protozoa AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Protozoa" ]]; then
        mkdir Analysis/LengthSeq/Protozoa
    fi

    if [[ -f "Analysis/LengthSeq/Protozoa/PROTOZOA_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Protozoa/PROTOZOA_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Protozoa/PROTOZOA_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Protozoa/Length_AA_PROTOZOA.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Protozoa AA Length analysis!\033[0m";
}
#
################## PLASTID ##################
#
function PLASTID_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Plastid DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Plastid" ]]; then
        mkdir Analysis/LengthSeq/Plastid
    fi

    if [[ -f "Analysis/LengthSeq/Plastid/PLASTID_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Plastid/PLASTID_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Plastid/PLASTID_DNA_ID.txt
            ./src/length_sequence.sh $file DNA UNZIPPED >> Analysis/LengthSeq/Plastid/Length_DNA_PLASTID.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Plastid DNA Length analysis!\033[0m";
}
function PLASTID_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Plastid AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Plastid" ]]; then
        mkdir Analysis/LengthSeq/Plastid
    fi

    if [[ -f "Analysis/LengthSeq/Plastid/PLASTID_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Plastid/PLASTID_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Plastid/PLASTID_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Plastid/Length_AA_PLASTID.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Plastid AA Length analysis!\033[0m";
}
#
# ################## MITOCHONDRIAL ##################
#
function MITOCHONDRIAL_DNA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Mitochondrial DNA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Mitochondrial" ]]; then
        mkdir Analysis/LengthSeq/Mitochondrial
    fi

    if [[ -f "Analysis/LengthSeq/Mitochondrial/MITOCHONDRIAL_DNA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Mitochondrial/MITOCHONDRIAL_DNA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Mitochondrial/MITOCHONDRIAL_DNA_ID.txt
            ./src/length_sequence.sh $file DNA UNZIPPED >> Analysis/LengthSeq/Mitochondrial/Length_DNA_MITOCHONDRIAL.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Mitochondrial DNA Length analysis!\033[0m";
}
function MITOCHONDRIAL_AA_LENGTH () {
    echo -e "\033[1;34m[RFSC]\033[0m Start Mitochondrial AA Length Analysis"

    if [[ ! -d "Analysis/LengthSeq" ]]; then
        echo "\033[1;34m[RFSC]\033[0m Folder Analysis/LengthSeq was not found! Aborting..."
        exit 0
    fi

    if [[ ! -d "Analysis/LengthSeq/Mitochondrial" ]]; then
        mkdir Analysis/LengthSeq/Mitochondrial
    fi

    if [[ -f "Analysis/LengthSeq/Mitochondrial/MITOCHONDRIAL_AA_ID.txt" ]]; then
        readarray -t already_processed < Analysis/LengthSeq/Mitochondrial/MITOCHONDRIAL_AA_ID.txt
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
            printf "$out_file\n" >> Analysis/LengthSeq/Mitochondrial/MITOCHONDRIAL_AA_ID.txt
            ./src/length_sequence.sh $file AA UNZIPPED >> Analysis/LengthSeq/Mitochondrial/Length_AA_MITOCHONDRIAL.txt
        fi
    done
    echo -e "\033[1;34m[RFSC]\033[1;32m Finished Mitochondrial AA Length analysis!\033[0m";
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_DNA_LENGTH;
  VIRUS_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Viral Length Sequences (DNA & AA) BD"
  ./src/length_sequence_csv.sh Virus VIRAL length_Viral.csv
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_DNA_LENGTH;
  BACTERIA_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Bacteria Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Bacteria BACTERIA length_Bacteria.csv
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_DNA_LENGTH;
  ARCHAEA_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Archaea Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Archaea ARCHAEA length_Archaea.csv
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_DNA_LENGTH;
  PROTOZOA_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Protozoa Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Protozoa PROTOZOA length_Protozoa.csv
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_DNA_LENGTH;
  FUNGI_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Fungi Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Fungi FUNGI length_Fungi.csv
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_DNA_LENGTH;
  PLANT_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plant Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Plant PLANT length_Plant.csv
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITOCHONDRIAL_DNA_LENGTH;
  MITOCHONDRIAL_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Mitochondrial Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Mitochondrial MITOCHONDRIAL length_Mitochondrial.csv
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_DNA_LENGTH;
  PLASTID_AA_LENGTH;
  echo -e "\033[1;34m[RFSC]\033[0m Start CSV generation for Plastid Length Sequences (DNA & AA) DB"
  ./src/length_sequence_csv.sh Plastid PLASTID length_Plastid.csv
  fi
#
# ==============================================================================