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
#
TRIMMING_FLAG=0;
TRIMMING_TYPE="";
TRIMMING_THREADS=0;
#
ASSEMBLY_FLAG=0;
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
# TRIMMING/FILTERING SEQUENCES
#
TRIMMING_SEQUENCE() {
	if [[ $TRIMMING_TYPE == "TT" ]]; then
		echo "Trimming with Trimmomatic"
		trimmomatic PE -threads 8 -phred33 SyntheticData/sample1.fq.gz SyntheticData/sample2.fq.gz GeneratedFiles/o_fw_pr.fq GeneratedFiles/o_fw_unpr.fq GeneratedFiles/o_rv_pr.fq GeneratedFiles/o_rv_unpr.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
	elif [[ $TRIMMING_TYPE == "FP" ]]; then
		echo "Trimming with FASTP"
		fastp -i SyntheticData/sample1.fq.gz -I SyntheticData/sample2.fq.gz -o GeneratedFiles/out1.fq.gz -O GeneratedFiles/out2.fq.gz
		mv fastp.html Outputs/
		mv fastp.json Outputs/
	else
		echo "Invalid Argument - $TRIMMING_TYPE!";
		echo "Use one of the follow:";
		echo "TT : To use the Trimmomatic Tool";
		echo "FP : To use the FASTP Tool";
		exit 0;
	fi
}
#
# ==================================================================
# DE-NOVO ASSEMBLY
#
SPADES_ASSEMBLY() {
	if [[ $TRIMMING_TYPE == "TT" ]]; then
		spades.py -t 16 --careful -o GeneratedFiles/out_spades_$1 -1 GeneratedFiles/o_fw_pr.fq -2 GeneratedFiles/o_rv_pr.fq -s GeneratedFiles/o_fw_unpr.fq -s GeneratedFiles/o_rv_unpr.fq
	else
		spades.py -t 16 --careful -o GeneratedFiles/out_spades_$1 -1 GeneratedFiles/out1.fq.gz -2 GeneratedFiles/out2.fq.gz 
	fi

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
		-trim|--filter)
			TRIMMING_FLAG=1;
			TRIMMING_TYPE="$2";
			TRIMMING_THREADS="$3";
			shift 3
		;;
		-rda|--run-de-novo)
			ASSEMBLY_FLAG=1;
			shift
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
	echo "                                                         "
	echo "   -trim, --filter <MODE> <THREADS>                      "
	echo "                  Filter Reads using Trimmomatic (TT)    "
	echo "                  or using FASTP (FP)                    "
	echo "                                                         "
	echo "   -rda, --run-de-novo                                   "
	echo "                  De-Novo Sequence Assembly              "
	echo "                                                         "
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
if [[ "$TRIMMING_FLAG" -eq "1" ]]; then
	echo "Start Trimming the Sequence!"
	TRIMMING_SEQUENCE "$TRIMMING_TYPE" "$TRIMMING_THREADS";
fi
#
# ===================================================================
#
if [[ "$ASSEMBLY_FLAG" -eq "1" ]]; then
	echo "Start De-Novo Assembly!"
	SPADES_ASSEMBLY;
fi
#
