#!/bin/bash

# Script used to analyse a given domain database outputing the 
# calculated predictions regarding the domain found
#
# Generates a file with the set of prediction in the /Predictions folder
#
# Run: ./Tests/deepTest.sh

SIZE=$1     # Percentage of the training dataset (ex: 0.8)

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
while [[ $# -gt 1 ]]
  do
  i="$2";
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
  echo "Usage: ./deepTest_GNB_TT.sh <Percentage> [options]                              "
  echo "                                                                                "
  echo "Script for training and assess the level of accuracy of the classifier          "
  echo "                                                                                "
  echo "   -h,   --help                   Show this help message,                       "
  echo "                                                                                "
  echo "   -vi,  --viral                  Run for viral DB,                             "
  echo "   -ba,  --bacteria               Run for bacteria DB,                          "
  echo "   -ar,  --archaea                Run for archaea DB,                           "
  echo "   -pr,  --protozoa               Run for protozoa DB,                          "
  echo "   -fu,  --fungi                  Run for fungi DB,                             "
  echo "   -pl,  --plant                  Run for plant DB,                             "
  echo "   -mi,  --mitochondrion          Run for mitochondrion DB,                     "
  echo "   -ps,  --plastid                Run for plastid DB.                           "
  echo "                                                                                "
  echo "Example: ./deepTest_GNB_TT.sh 0.8 --viral --fungi --plant                       "
  echo "                                                                                "
  exit 1
  fi
#
# ==============================================================================
#
################## VIRUS ##################
#
function VIRUS_TEST () {
    virus_files_NM=`wc -l Analysis/GeCo/Virus/geco3_Viral.csv | awk '{ print $1 }'`

    training_seq=`echo $virus_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$virus_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Viral...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)

    done <Analysis/GCcontent/Virus/gc_content_Viral.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)

    done <Analysis/GeCo/Virus/geco3_Viral.csv

    while read line; 
    do
        AA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)

    done <Analysis/AC/Virus/ac2_Viral.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Virus/length_Viral.csv

    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Viral.txt
        fi
    done
}
#
################## BACTERIA ##################
#
function BACTERIA_TEST () {
    bacteria_files_NM=`wc -l Analysis/GeCo/Bacteria/geco3_Bacteria.csv | awk '{ print $1 }'`

    training_seq=`echo $bacteria_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$bacteria_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Bacteria...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Bacteria/gc_content_Bacteria.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Bacteria/geco3_Bacteria.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Bacteria/ac2_Bacteria.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Bacteria/length_Bacteria.csv

    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Bacteria.txt
        fi
    done
}
#
################## ARCHAEA ##################
#
function ARCHAEA_TEST () {
    archaea_files_NM=`wc -l Analysis/GeCo/Archaea/geco3_Archaea.csv | awk '{ print $1 }'`

    training_seq=`echo $archaea_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$archaea_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Archaea...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Archaea/gc_content_Archaea.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Archaea/geco3_Archaea.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Archaea/ac2_Archaea.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Archaea/length_Archaea.csv


    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Archaea.txt
        fi
    done
}
#
################## FUNGI ##################
#
function FUNGI_TEST () {
    fungi_files_NM=`wc -l Analysis/GeCo/Fungi/geco3_Fungi.csv | awk '{ print $1 }'`

    training_seq=`echo $fungi_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$fungi_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Fungi...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Fungi/gc_content_Fungi.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Fungi/geco3_Fungi.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Fungi/ac2_Fungi.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Fungi/length_Fungi.csv

    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Fungi.txt
        fi
    done
}
#
################## PLANT ##################
#
function PLANT_TEST () {
    plant_files_NM=`wc -l Analysis/GeCo/Plant/geco3_Plant.csv | awk '{ print $1 }'`

    training_seq=`echo $plant_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$plant_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Plant...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Plant/gc_content_Plant.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Plant/geco3_Plant.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Plant/ac2_Plant.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Plant/length_Plant.csv


    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Plant.txt
        fi
    done
}
#
################## PROTOZOA ##################
#
function PROTOZOA_TEST () {
    protozoa_files_NM=`wc -l Analysis/GeCo/Protozoa/geco3_Protozoa.csv | awk '{ print $1 }'`

    training_seq=`echo $protozoa_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$protozoa_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Protozoa...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Protozoa/geco3_Protozoa.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Protozoa/ac2_Protozoa.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Protozoa/length_Protozoa.csv


    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Protozoa.txt
        fi
    done
}
#
################## MITOCHONDRIAL ##################
#
function MITO_TEST () {
    mito_files_NM=`wc -l Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv | awk '{ print $1 }'`

    training_seq=`echo $mito_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$mito_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Mitochondrial...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Mitochondrial/gc_content_Mitochondrial.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Mitochondrial/ac2_Mitochondrial.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Mitochondrial/length_Mitochondrial.csv


    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Mitochondrial.txt
        fi
    done
}
#
################## PLASTID ##################
#
function PLASTID_TEST () {
    plastid_files_NM=`wc -l Analysis/GeCo/Plastid/geco3_Plastid.csv | awk '{ print $1 }'`

    training_seq=`echo $plastid_files_NM*$SIZE -1 | bc`
    training_seq=${training_seq%.*}

    START=1;
    STOP=$plastid_files_NM;

    echo -e "\033[1;34m[RFSC]\033[0m Analysing Plastid...";

    while read line; 
    do
        GC_VALUE+=("0"`echo $line | cut -d ' ' -f4`)
        
    done <Analysis/GCcontent/Plastid/gc_content_Plastid.csv

    while read line; 
    do
        DNA_VALUE_1+=("0"`echo $line | cut -d ' ' -f7`)
        DNA_VALUE_3+=("0"`echo $line | cut -d ' ' -f9`)
        DNA_VALUE_7+=("0"`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/GeCo/Plastid/geco3_Plastid.csv

    while read line; 
    do
        AA_VALUE_1+=(`echo $line | cut -d ' ' -f7`)
        AA_VALUE_3+=(`echo $line | cut -d ' ' -f9`)
        AA_VALUE_7+=(`echo $line | cut -d ' ' -f11`)
        
    done <Analysis/AC/Plastid/ac2_Plastid.csv

    while read line; 
    do
        DNA_LENGTH+=(`echo $line | cut -d ' ' -f2`)
        AA_LENGTH+=(`echo $line | cut -d ' ' -f3`)
        
    done <Analysis/LengthSeq/Plastid/length_Plastid.csv


    for seq in $(eval echo "{$START..$STOP}")
    do
        if [[ "$seq" -gt "${training_seq}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/naiveBayes.py 8 $SIZE 0 0 ${DNA_VALUE_1[seq]} ${AA_VALUE_1[seq]} ${DNA_VALUE_3[seq]} ${AA_VALUE_3[seq]} ${DNA_VALUE_7[seq]} ${AA_VALUE_7[seq]} ${GC_VALUE[seq]} ${DNA_LENGTH[seq]} ${AA_LENGTH[seq]} 1111 >> Tests/Predictions/GNB/Test_Train_Database/Prediction_Plastid.txt
        fi
    done
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRUS_TEST;
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_TEST;
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_TEST;
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_TEST;
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_TEST;
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_TEST;
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITO_TEST;
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_TEST;
  fi
#
# ==============================================================================
#