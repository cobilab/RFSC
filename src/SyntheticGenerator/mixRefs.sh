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
        echo -e "\033[1;34m[NCRS]\033[0m Start File Concat..."
else
	echo -e "\033[1;34m[NCRS] \033[1;31m Argument Error - Insert Reference(s) File(s) \033[0m"
        echo -e "\033[1;34m[NCRS] \033[0;33m Usage: bash mixRefs.sh [FILE1] [FILE2] : [FILE9] \033[0m"
	exit 0
fi

NUM_FILES=$( ls src/SyntheticGenerator/Inputs/ | wc -l )
NUM_FILES=$((NUM_FILES+1))

zcat $1 $2 $3 $4 $5 $6 $7 $8 $9 > "src/SyntheticGenerator/Inputs/mix${NUM_FILES}.fa"
echo -e "\033[1;34m[NCRS] \033[0;32m Success! File mix${NUM_FILES}.fa has been created. \033[0m"