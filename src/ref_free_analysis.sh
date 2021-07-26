#!/bin/bash

METHOD=$1;          # GNB | KNN
#
# KNN Parameters
K=$2;               # K-Neighbors
#
# GNB Parameters
N_DOMAINS=$3;       # Number of Domains
PREDICTOR_CODE=$4;  # Default: All predictors: 1111
#
if [ -f Input_Data/ReferenceFree/GeneratedPredictorValues/*.csv ]; then
    NUMBER_INPUTS=`wc -l Input_Data/ReferenceFree/GeneratedPredictorValues/Input_Predictors.csv | awk '{ print $1 }'`
    NUMBER_INPUTS=$(( NUMBER_INPUTS - 1 ))
else
    echo -e "\033[1;34m[RFSC]\033[0m Analyzes of predictors of input sequences not found"
    echo -e "\033[1;34m[RFSC]\033[0m Run the ./src/get_input_predictors.sh script"
    exit 0
fi

if [[ ! -d "Results/ReferenceFree" ]]; then
    mkdir Results/ReferenceFree
    echo -e "\033[1;34m[RFSC]\033[0m Alert: The results will be stored at Results/ReferenceFree"
fi

if [ "$METHOD" = "GNB" ]; then
    echo -e "\033[1;34m[RFSC]\033[0m Starting Gaussian Naive Bayes Classification"

    while read line; 
    do
        IDENTIFIER+=(`echo $line | cut -d ',' -f1`)
        GECO3_1+=(`echo $line | cut -d ',' -f2`)
        AC_1+=(`echo $line | cut -d ',' -f3`)
        GECO3_3+=(`echo $line | cut -d ',' -f4`)
        AC_3+=(`echo $line | cut -d ',' -f5`)
        GECO3_7+=(`echo $line | cut -d ',' -f6`)
        AC_7+=(`echo $line | cut -d ',' -f7`)
        GC+=(`echo $line | cut -d ',' -f8`)
        DNA_LEN+=(`echo $line | cut -d ',' -f9`)
        AA_LEN+=(`echo $line | cut -d ',' -f10`)
    done <Input_Data/ReferenceFree/GeneratedPredictorValues/Input_Predictors.csv

    START=1;
    STOP=$NUMBER_INPUTS;

    for seq in $(eval echo "{$START..$STOP}")
    do
        echo -e "\033[1;34m[RFSC]\033[0m Analysing: ${IDENTIFIER[seq]}";
        PREDICTION=`python3 src/naiveBayes.py $N_DOMAINS 0 1 1 ${GECO3_1[seq]} ${AC_1[seq]} ${GECO3_3[seq]} ${AC_3[seq]} ${GECO3_7[seq]} ${AC_7[seq]} ${GC[seq]} ${DNA_LEN[seq]} ${AA_LEN[seq]} $PREDICTOR_CODE`
        echo "${IDENTIFIER[seq]} -> $PREDICTION" >> Results/ReferenceFree/GNB_Predictions.txt
        echo -e "\033[1;34m[RFSC]\033[0m Using Gaussian Naive Bayes classification ${IDENTIFIER[seq]} has been classified as $PREDICTION";
    done

elif [ "$METHOD" = "KNN" ]; then
    echo -e "\033[1;34m[RFSC]\033[0m Starting K-Nearest Neighbor Classification"

    while read line; 
    do
        IDENTIFIER+=(`echo $line | cut -d ',' -f1`)
        GECO3_3+=(`echo $line | cut -d ',' -f4`)
        AC_3+=(`echo $line | cut -d ',' -f5`)
        GC+=(`echo $line | cut -d ',' -f8`)
        DNA_LEN+=(`echo $line | cut -d ',' -f9`)
        AA_LEN+=(`echo $line | cut -d ',' -f10`)
    done <Input_Data/ReferenceFree/GeneratedPredictorValues/Input_Predictors.csv

    START=1;
    STOP=$NUMBER_INPUTS;

    for seq in $(eval echo "{$START..$STOP}")
    do
        echo -e "\033[1;34m[RFSC]\033[0m Analysing: ${IDENTIFIER[seq]}";
        PREDICTION=`python3 src/KNN.py Test $K 0 1 ${GECO3_3[seq]} ${AC_3[seq]} ${GC[seq]} ${DNA_LEN[seq]} ${AA_LEN[seq]} NULL`
        echo "${IDENTIFIER[seq]} -> $PREDICTION" >> Results/ReferenceFree/KNN_Predictions.txt
        echo -e "\033[1;34m[RFSC]\033[0m Using K-Nearest Neighbor classification ${IDENTIFIER[seq]} has been classified as $PREDICTION";
    done

else
    echo -e "\033[1;34m[RFSC]\033[0m Classification Method not recognised! \n\033[1;34m[RFSC]\033[0m Please choose a method between GNB and KNN"
fi