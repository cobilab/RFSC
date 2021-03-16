#!/bin/bash

# ===================================================== #
# =                                                   = #
# =                  References Mixer                 = #
# =                                                   = #
# ===================================================== #
# =                                                   = #
# =   Usage: sh mixRefs.sh [FILE1] [FILE2] : [FILE9]  = #
# =                                                   = #
# ===================================================== #

# VARIABLES
NUM_FILES=0;
FILE1=$1;
FILE2=$2;
FILE3=$3;
FILE4=$4;
FILE5=$5;
FILE6=$6;
FILE7=$7;
FILE8=$8;
FILE9=$9;

if [ -n "$FILE1" ]; then
        echo "Start File Concat..."
else
	echo "Argument Error - Insert Reference(s) File(s)"
        echo "Usage: bash mixRefs.sh [FILE1] [FILE2] : [FILE9]"
	exit 0
fi

NUM_FILES=$( ls SyntheticGenerator/Inputs/ | wc -l )
NUM_FILES=$((NUM_FILES+1))

zcat $1 $2 $3 $4 $5 $6 $7 $8 $9 > "SyntheticGenerator/Inputs/mix${NUM_FILES}.fa"
echo "Success! File mix${NUM_FILES}.fa has been created."

#echo "${NUM_FILES}"
