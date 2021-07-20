#!/bin/bash

# Script used to analyse a given domain database outputing the 
# calculated predictions regarding the domain found
#
# Generates a file with the set of prediction in the /Predictions/KNN folder
#
# Run: ./Tests/deepTest_KNN.sh

blocks_perc=(0 0.2 0.4 0.6 0.8 1) # Five-Fold Cross Validation (5 blocks)
num_blocks=$(( ${#blocks_perc[@]} - 2 ))
#
K=2
#
TEST=0;
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
COUNT_SAMPLE_VIRUS=0;
COUNT_SAMPLE_BACT=0;
COUNT_SAMPLE_ARCH=0;
COUNT_SAMPLE_FUNGI=0;
COUNT_SAMPLE_PLANT=0;
COUNT_SAMPLE_PROTO=0;
COUNT_SAMPLE_MITO=0;
COUNT_SAMPLE_PLASTID=0;
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
    -t|--test)                  TEST=1;                 HELP=0; shift  ;;
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
  echo "                                                                                         "
  echo "Usage: ./deepTest_KNN.sh [options]                                                       "
  echo "                                                                                         "
  echo "Script for training and assess the level of accuracy of the classifier                   "
  echo "                                                                                         "
  echo "   -h,   --help                   Show this help message,                                "
  echo "                                                                                         "
  echo "   -t,   --test                   Run Test Database against Train Database               "
  echo "                                                                                         "
  echo "   -vi,  --viral                  Run for Train viral DB, with Cross-Validation,         "
  echo "   -ba,  --bacteria               Run for Train bacteria DB, with Cross-Validation,      "
  echo "   -ar,  --archaea                Run for Train archaea DB, with Cross-Validation,       "
  echo "   -pr,  --protozoa               Run for Train protozoa DB, with Cross-Validation,      "
  echo "   -fu,  --fungi                  Run for Train fungi DB, with Cross-Validation,         "
  echo "   -pl,  --plant                  Run for Train plant DB, with Cross-Validation,         "
  echo "   -mi,  --mitochondrion          Run for Train mitochondrion DB, with Cross-Validation, "
  echo "   -ps,  --plastid                Run for Train plastid DB, with Cross-Validation.       "
  echo "                                                                                         "
  echo "Example: ./deepTest_KNN.sh --test                                                        "
  echo "                                                                                         "
  exit 1
  fi
#
# ==============================================================================
#
function TEST_DATABASE () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing Test Database...";
    readarray -t trainDatabase < Analysis/KNN/Test.csv

    for sample in "${trainDatabase[@]:1}"
    do
        #echo "${sample}"
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA=(`echo ${sample} | cut -d ',' -f2`)
        AA=(`echo ${sample} | cut -d ',' -f3`)
        GC=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA=(`echo ${sample} | cut -d ',' -f5`)
        L_AA=(`echo ${sample} | cut -d ',' -f6`)

        echo -e "\033[1;34m[RFSC]\033[0m Run: $DOMAIN";
        PREDICTION=`python3 src/KNN.py Test $K 0 1 $DNA $AA $GC $L_DNA $L_AA NULL`
        echo "$DOMAIN,$PREDICTION" >> Tests/Predictions/KNN/TestDatabase.txt
    done
}
#
################## CROSS VALIDATION ##################
function VIRAL_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Viral Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=1;
    STOP=$COUNT_SAMPLE_VIRUS;

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_VIRUS*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_VIRUS*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Viral >> Tests/Predictions/KNN/Prediction_Viral_CV.txt
            fi
        done
    done
}
#
function BACTERIA_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Bacteria Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+1));
    STOP=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_BACT*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_BACT*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Bacteria >> Tests/Predictions/KNN/Prediction_Bacteria_CV.txt
            fi
        done
    done
}
#
function ARCHAEA_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Archaea Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+1))
    STOP=$(($COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_ARCH*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_ARCH*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Archaea >> Tests/Predictions/KNN/Prediction_Archaea_CV.txt
            fi
        done
    done
}
#
function FUNGI_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Fungi Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+1))
    STOP=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_FUNGI*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_FUNGI*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Fungi >> Tests/Predictions/KNN/Prediction_Fungi_CV.txt
            fi
        done
    done
}
#
function PLANT_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Plant Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+1))
    STOP=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_PLANT*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_PLANT*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Plant >> Tests/Predictions/KNN/Prediction_Plant_CV.txt
            fi
        done
    done
}
#
function PROTOZOA_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Protozoa Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT+1))
    STOP=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT+$COUNT_SAMPLE_PROTO));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_PROTO*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_PROTO*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Protozoa >> Tests/Predictions/KNN/Prediction_Protozoa_CV.txt
            fi
        done
    done
}
#
function MITO_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Mitochondrial Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT+$COUNT_SAMPLE_PROTO+1))
    STOP=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT+$COUNT_SAMPLE_PROTO+$COUNT_SAMPLE_MITO));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_MITO*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_MITO*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Mitochondrial >> Tests/Predictions/KNN/Prediction_Mitochondrial_CV.txt
            fi
        done
    done
}
#
function PLASTID_CV () {
    echo -e "\033[1;34m[RFSC]\033[0m Analysing KNN Plastid Cross-Validation...";

    readarray -t trainDatabase < Analysis/KNN/Train.csv

    for sample in "${trainDatabase[@]:1}"
    do
        DOMAIN=(`echo ${sample} | cut -d ',' -f1`)
        DNA+=(`echo ${sample} | cut -d ',' -f2`)
        AA+=(`echo ${sample} | cut -d ',' -f3`)
        GC+=(`echo ${sample} | cut -d ',' -f4`)
        L_DNA+=(`echo ${sample} | cut -d ',' -f5`)
        L_AA+=(`echo ${sample} | cut -d ',' -f6`)

        if [[ $DOMAIN =~ "Viral" ]]; then
            (( COUNT_SAMPLE_VIRUS++ ))
        elif [[ $DOMAIN =~ "Bacteria" ]]; then
            (( COUNT_SAMPLE_BACT++ ))
        elif [[ $DOMAIN =~ "Archaea" ]]; then
            (( COUNT_SAMPLE_ARCH++ ))
        elif [[ $DOMAIN =~ "Fungi" ]]; then
            (( COUNT_SAMPLE_FUNGI++ ))
        elif [[ $DOMAIN =~ "Plant" ]]; then
            (( COUNT_SAMPLE_PLANT++ ))
        elif [[ $DOMAIN =~ "Protozoa" ]]; then
            (( COUNT_SAMPLE_PROTO++ ))
        elif [[ $DOMAIN =~ "Mitochondrial" ]]; then
            (( COUNT_SAMPLE_MITO++ ))
        elif [[ $DOMAIN =~ "Plastid" ]]; then
            (( COUNT_SAMPLE_PLASTID++ ))
        fi
    done

    START=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT+$COUNT_SAMPLE_PROTO+$COUNT_SAMPLE_MITO+1))
    STOP=$(($COUNT_SAMPLE_VIRUS+$COUNT_SAMPLE_BACT+$COUNT_SAMPLE_ARCH+$COUNT_SAMPLE_FUNGI+$COUNT_SAMPLE_PLANT+$COUNT_SAMPLE_PROTO+$COUNT_SAMPLE_MITO+$COUNT_SAMPLE_PLASTID));

    for i in $(eval echo "{0..$num_blocks}")
    do
        training_block_min=`echo $COUNT_SAMPLE_PLASTID*${blocks_perc[i]} | bc`
        training_block_min=${training_block_min%.*}
        training_block_min=$(($START+$training_block_min))
        training_block_max=`echo $COUNT_SAMPLE_PLASTID*${blocks_perc[i+1]} | bc`
        training_block_max=${training_block_max%.*}
        training_block_max=$(($START+$training_block_max))

        echo -e "\033[1;34m[RFSC]\033[0m Testing sequences from: $training_block_min - $training_block_max:"

        for seq in $(eval echo "{$START..$STOP}")
        do
            if [[ "$seq" -ge "${training_block_min}" ]] && [[ "$seq" -le "${training_block_max}" ]]; then
            echo -e "\033[1;34m[RFSC]\033[0m Run: $seq";
            python3 src/KNN.py CV $K ${blocks_perc[i]} ${blocks_perc[i+1]} ${DNA[seq]} ${AA[seq]} ${GC[seq]} ${L_DNA[seq]} ${L_AA[seq]} Plastid >> Tests/Predictions/KNN/Prediction_Plastid_CV.txt
            fi
        done
    done
}
#
# ==============================================================================
#
if [[ "$TEST" -eq "1" ]];
  then
  TEST_DATABASE;
  fi
#
if [[ "$VIRAL" -eq "1" ]];
  then
  VIRAL_CV;
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BACTERIA_CV;
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  ARCHAEA_CV;
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  PROTOZOA_CV;
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  FUNGI_CV;
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  PLANT_CV;
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  MITO_CV;
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  PLASTID_CV;
  fi
#
# ==============================================================================
#