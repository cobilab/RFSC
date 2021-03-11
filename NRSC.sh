#!/bin/bash
#
#####################################################################
# ================================================================= #
# =                                                               = #
# =                            N R S C                            = #
# =								  = #
# =         A Non-Referencial Sequence Classification Tool        = # 
# =            for DNA sequences in metagenomic samples.          = #
# =                                                               = #
# ================================================================= #
#####################################################################
#
SHOW_HELP=0;
SHOW_VERSION=0;
#
GEN_SYNTHETIC=0;
REF_FILE1="";
REF_FILE2="";
REF_FILE3="";
REF_FILE4="";
LAST_MIX=0;
#
CRYFA_FLAG=0;
TRIMMOMATIC_FLAG=0;
FASTP_FLAG=0;
FALCON_FLAG=0;
BLASTN_FLAG=0;
#
# ==================================================================
# CURRENT VIRUSES OR VIRUSES GROUPS ACCEPTED TO BE SEARCHED
#
declare -a VIRUSES=("B19" "HBV");
#
# ==================================================================
# GENERATE SYNTHETIC SEQUENCE
#
GENERATE_SYNTHETIC () {
	./SyntheticGenerator/mixRefs.sh "$REF_FILE1" "$REF_FILE2" "$REF_FILE3";
	LAST_MIX=$( ls SyntheticGenerator/Inputs/ | wc -l )
	./SyntheticGenerator/syntheticGen.sh SyntheticGenerator/Inputs/"mix${LAST_MIX}.fa";
}
#
# ==================================================================
# OPTIONS
#
if [ "$#" -eq 0 ]; then
	SHOW_HELP=1;
fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
	do
	i="$1";
	case $i in
		-h|--help|?)
			SHOW_HELP=1;
			shift
		;;
		-v|-V|--version)
			SHOW_VERSION=1;
			shift
		;;
		-synt|--synthetic)
			GEN_SYNTHETIC=1;
			REF_FILE1="$2";
			REF_FILE2="$3";
			REF_FILE3="$4";
			shift 4
		;;
		-*) # Unknown option
		echo "Invalid arg ($1)!";
		echo "For help, try: bash NRSC.sh -h"
		exit 1;
		;;
	esac
done
#
set -- "${POSITIONAL[@]}" # Restore positional parameters
#
# ======================================================================
# HELP MENU
#
if [ "$SHOW_HELP" -eq "1" ]; then

	echo "                                                         "
	echo "                 _   _ ___   ____   ____                 "
	echo "                | \ | |  _ \/ ___| / ___|                "
	echo "                |  \| | |_) \___ \| |                    "
	echo "                | |\  |  _ < ___) | |___                 "
	echo "                |_| \_|_| \_\____/ \____|                "
	echo "                                                         "
	echo "                     P I P E L I N E                     "
	echo "                                                         "
	echo "      A Non-Referencial Sequence Classification Tool     "
	echo "        for DNA sequences in metagenomic samples.        "
	echo "                                                         "
	echo "                Usage: ./NRSC.sh [options]               "
	echo "                                                         "
	echo "   -h,  --help	Show this help message and exit,       "
	echo "   -v,  --version Show the version and some information  "
	echo "                                                         "
	echo "   -synt, --synthetic [FILE1] : [FILE3]                  "
	echo "                  Generate a synthetical sequence using  "
	echo "                  3 reference files for testing purposes "
	exit 1;
fi
#
# ======================================================================
# VERSION
#
if [ "$SHOW_VERSION" -eq "1" ]; then
	
	echo "                                                         "
	echo "                          NRSC                           "
	echo "                                                         "
	echo "                      Version: 1.0                       "
	echo "                                                         "
	echo "                       IEETA/DETI                        "
	echo "             University of Aveiro, Portugal.             "
	echo "                                                         "
	exit 0;
fi	
#
# ===================================================================
#
if [[ "$GEN_SYNTHETIC" -eq "1" ]]; then
	echo "Start Synthetic Sequence Generation!"
	GENERATE_SYNTHETIC "$REF_FILE1" "$REF_FILE2" "$REF_FILE3";
fi
#
# ===================================================================
#
