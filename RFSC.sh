#!/bin/bash
#
#####################################################################
# ================================================================= #
# =                                                               = #
# =                            R F S C                            = #
# =                                                               = #
# =         A Reference-Free Sequence Classification Tool         = # 
# =            for DNA sequences in metagenomic samples.          = #
# =                                                               = #
# ================================================================= #
#####################################################################
#
SHOW_HELP=0;
SHOW_VERSION=0;
INSTALL=0;
#
BUILD_DB_VIRUS=0;
#
GEN_ADAPTERS=0;
#
THREADS_AVAILABLE=`cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l`;
THREADS=0;
GET_THREADS="";
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
TRIMMING_MODE=""; # SE or PE
#
ASSEMBLY_FLAG=0;
#
FALCON_FLAG=0;
FALCON_MODE="";
BLASTN_FLAG=0;
#
RUN_DECRYPT=0;
RUN_ENCRYPT=0;
#
# ==================================================================
# CURRENT VIRUSES OR VIRUSES GROUPS ACCEPTED TO BE SEARCHED
#
declare -a VIRUSES=("B19" "HBV");
#
# ==================================================================
# VERIFICATION FUNCTIONS
#
CHECK_ADAPTERS() {
	if [ ! -f Input_Data/adapters.fa ]; then
		echo -e "\033[1;33m[RFSC] ERROR: adapter sequences (adapter.fa) not found! \033[0m"
		echo -e "\033[1;34m[RFSC] \033[0;33m ./RFSC.sh --gen-adapters \033[0m : To generate the adapter sequences ...";
	fi
}
#
# ==================================================================
# GENERATE SYNTHETIC SEQUENCE
#
GENERATE_SYNTHETIC () { 
	./src/SyntheticGenerator/mixRefs.sh "$REF_FILE1" "$REF_FILE2" "$REF_FILE3";
	LAST_MIX=$( ls src/SyntheticGenerator/Inputs/ | wc -l )
	./src/SyntheticGenerator/syntheticGen.sh src/SyntheticGenerator/Inputs/"mix${LAST_MIX}.fa";
}
#
# ==================================================================
# TRIMMING/FILTERING SEQUENCES
#
TRIMMING_SEQUENCE() {
	#
	# Fetch the input files
	i=0;
	for file in Input_Data/*.fq.gz
	do
		input_file[i]="$file"
		(( i++ ))
	done
	#
	if [[ $TRIMMING_TYPE == "TT" ]]; then
		echo -e "\033[1;34m[RFSC]\033[0m Trimming using Trimmomatic";
		TRIMMING_THREADS=$THREADS_AVAILABLE;
		echo -e "\033[1;34m[RFSC]\033[0m Currently using $TRIMMING_THREADS available threads!";
		#
		CHECK_ADAPTERS;
		cp Input_Data/adapters.fa adapters.fa
		#
		if [[ $TRIMMING_MODE == "PE" ]]; then
			trimmomatic $TRIMMING_MODE -threads $TRIMMING_THREADS -phred33 ${input_file[0]} ${input_file[1]} GeneratedFiles/o_fw_pr.fq GeneratedFiles/o_fw_unpr.fq GeneratedFiles/o_rv_pr.fq GeneratedFiles/o_rv_unpr.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
		elif [[ $TRIMMING_MODE == "SE" ]]; then
			trimmomatic $TRIMMING_MODE -threads $TRIMMING_THREADS -phred33 ${input_file[0]} GeneratedFiles/out.fq ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:25
		else
			echo -e "\033[1;34m[RFSC] \033[1;31m Invalid Argument - $TRIMMING_MODE! \033[0m";
			echo -e "\033[1;34m[RFSC]\033[0m Use one of the follow:";
			echo -e "\033[1;34m[RFSC] \033[0;33m PE \033[0m : To use paired end reads";
			echo -e "\033[1;34m[RFSC] \033[0;33m SE \033[0m : To use single end reads";
			exit 0;
		fi
		#
		rm adapters.fa
		#
	elif [[ $TRIMMING_TYPE == "FP" ]]; then
		echo -e "\033[1;34m[RFSC]\033[0m Trimming with FASTP";
		#
		if [[ $TRIMMING_MODE == "PE" ]]; then
			fastp -i ${input_file[0]} -I ${input_file[1]} -o GeneratedFiles/out1.fq.gz -O GeneratedFiles/out2.fq.gz
		elif [[ $TRIMMING_MODE == "SE" ]]; then
			fastp -i ${input_file[0]} -o GeneratedFiles/out.fq.gz
		else
			echo -e "\033[1;34m[RFSC] \033[1;31m Invalid Argument - $TRIMMING_MODE! \033[0m";
			echo -e "\033[1;34m[RFSC]\033[0m Use one of the follow:";
			echo -e "\033[1;34m[RFSC] \033[0;33m PE \033[0m : To use paired end reads";
			echo -e "\033[1;34m[RFSC] \033[0;33m SE \033[0m : To use single end reads";
			exit 0;
		fi
		#
		mv fastp.html Outputs/
		mv fastp.json Outputs/
	else
		echo -e "\033[1;34m[RFSC] \033[1;31m Invalid Argument - $TRIMMING_TYPE! \033[0m";
		echo -e "\033[1;34m[RFSC]\033[0m Use one of the follow:";
		echo -e "\033[1;34m[RFSC] \033[0;33m TT \033[0m : To use the Trimmomatic Tool";
		echo -e "\033[1;34m[RFSC] \033[0;33m FP \033[0m : To use the FASTP Tool";
		exit 0;
	fi
}
#
# ==================================================================
# DE-NOVO ASSEMBLY
#
SPADES_ASSEMBLY() {
	if [[ $TRIMMING_TYPE == "TT" ]]; then
		if [[ $TRIMMING_MODE == "PE" ]]; then
			spades.py -t 16 --careful -o GeneratedFiles/out_spades_$1 -1 GeneratedFiles/o_fw_pr.fq -2 GeneratedFiles/o_rv_pr.fq -s GeneratedFiles/o_fw_unpr.fq -s GeneratedFiles/o_rv_unpr.fq
		else
			spades.py -t 16 --careful -o GeneratedFiles/out_spades_$1 -s GeneratedFiles/out.fq
		fi
	else
		if [[ $TRIMMING_MODE == "PE" ]]; then
			spades.py -t 16 --careful -o GeneratedFiles/out_spades_$1 -1 GeneratedFiles/out1.fq.gz -2 GeneratedFiles/out2.fq.gz
		else
			spades.py -t 16 --careful -o GeneratedFiles/out_spades_$1 -s GeneratedFiles/out.fq.gz
		fi
	fi

}
#
# ==================================================================
# FALCON ANALYSIS - SCAFFOLDS
#
FALCON_SO_MODE(){
	FALCON -n $THREADS_AVAILABLE -v -F -x Outputs/falcon_SO_results.txt GeneratedFiles/out_spades_/scaffolds.fasta References/NCBI-Virus/VDB.fa
        echo -e "\033[1;34m[RFSC]\033[0m Outputs/falcon_SO_results.txt file was successfully been generated!"
	#
	readarray -t results <Outputs/falcon_SO_results.txt
	#
	NUM_RES=${#results[@]}
	for (( i=0; i<$NUM_RES; i++ ));
	do
		VIRUS=`echo "${results[i]}"|awk 'NF>1{print $NF}'`
		PER=`echo "${results[i]}"|awk '{print $3}'`
		BOOL=`echo "$PER > 70.000" | bc`
		#
		if [[ $BOOL -eq "1" ]]; then
				printf "$VIRUS\n" >> Results/ref_result.txt
		fi
	done
}
#
# ==================================================================
# FALCON ANALYSIS - EACH READS
#
FALCON_RM_MODE(){
	echo -e "\033[1;34m[RFSC]\033[0m Start Breaking File into Reads"
        mkdir GeneratedFiles/out_spades_/Reads
        awk '/>/{filename="GeneratedFiles/out_spades_/Reads/"NR".fasta"}; {print >filename}' GeneratedFiles/out_spades_/scaffolds.fasta
        # 
        reads=0
        while read line # Get all the Read Files Names
        do
        	array[ $reads ]="$line"
                (( reads++ ))
        done < <(ls -ls GeneratedFiles/out_spades_/Reads)
        #
        mkdir Outputs/FalconReads
        len=${#array[@]}
        i=1
		while [[ $i != $len ]]
		do
			if [[ $(( ($len - $i) % 2)) == 0 ]]; then
				R1=`echo "${array[i]}"|awk 'NF>1{print $NF}'`
				R2=`echo "${array[i+1]}"|awk 'NF>1{print $NF}'`
				echo -e "\033[1;34m[RFSC]\033[0m Analysing Read $R1 & $R2 with $(($THREADS_AVAILABLE/2)) threads each!"
				FALCON -n $(($THREADS_AVAILABLE/2)) -v -F -x Outputs/FalconReads/falcon_RM_"${R1}"_results.txt GeneratedFiles/out_spades_/Reads/$R1 References/NCBI-Virus/VDB.fa | FALCON -n $(($THREADS_AVAILABLE/2)) -v -F -x Outputs/FalconReads/falcon_RM_"${R2}"_results.txt GeneratedFiles/out_spades_/Reads/$R2 References/NCBI-Virus/VDB.fa
				(( i+=2 ))
			else
				R=`echo "${array[i]}"|awk 'NF>1{print $NF}'`
				echo -e "\033[1;34m[RFSC]\033[0m Analysing Read $R1 with $(($THREADS_AVAILABLE)) threads!"
				FALCON -n $(($THREADS_AVAILABLE)) -v -F -x Outputs/FalconReads/falcon_RM_"${R}"_results.txt GeneratedFiles/out_spades_/Reads/$R References/NCBI-Virus/VDB.fa
				(( i++ ))
			fi 
		done
}
#
# ==================================================================
# FALCON ANALYSIS
#
FALCON_ANALYSIS() {
	if [[ $FALCON_MODE == "SO" ]]; then
		FALCON_SO_MODE;
	elif [[ $FALCON_MODE == "RM" ]]; then
		FALCON_SO_MODE;
		FALCON_RM_MODE;
	else
		echo -e "\033[1;34m[RFSC] \033[1;31m Invalid Argument - $FALCON_MODE! \033[0m";
		echo -e "\033[1;34m[RFSC]\033[0m Use one of the follow:";
		echo -e "\033[1;34m[RFSC] \033[0;33m SO \033[0m : To analyse only the scaffold";
		echo -e "\033[1;34m[RFSC] \033[0;33m RM \033[0m : To use redundancy in the analyse (Run each Read)";
		exit 0;
	fi

}
#
# ==================================================================
# ENCRYPTION
#
ENCRYPT_DATA() {
	echo -e "\033[1;34m[RFSC]\033[0;32m Please insert the password to encrypt the files contained in /Data_Security/Encrypt_Input: \033[0m";
	read -s password
	echo "$password" > key.txt
	#
	for file in Data_Security/Encrypt_Input/*
	do
		echo -e "\033[1;34m[RFSC]\033[0;32m Encrypting $file ... \033[0m";
		out_file=$(basename $file);
		cryfa -k key.txt $file > Data_Security/Encrypted_Data/$out_file.enc
	done
	rm -f key.txt
	echo -e "\033[1;34m[RFSC]\033[0;32m $file has been encrypted! \033[0m";
}
#
# ==================================================================
# DECRYPTION
#
DECRYPT_DATA() {
	for file in Data_Security/Encrypted_Data/*
	do
		echo -e "\033[1;34m[RFSC]\033[0;32m Decrypting $file ... \033[0m";
		out_file=$(basename $file);
		echo -e "\033[1;34m[RFSC]\033[0;32m Please insert the password to decrypt the file $file: \033[0m";
		read -s password
		echo "$password" > key.txt
		cryfa -k key.txt -d $file > Data_Security/Decrypted_Data/$out_file.dec;
		echo -e "\033[1;34m[RFSC]\033[0;32m $file has been decrypted! \033[0m";
		rm -f key.txt;
	done
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
		-i|--install)
			INSTALL=1;
			shift
		;;
		-t|--threads)
			THREADS=1;
			GET_THREADS="$2";
			shift 2
		;;
		-bref|--build-ref-virus)
			BUILD_DB_VIRUS=1;
			shift
		;;
		-gad|--gen-adapters)
			GEN_ADAPTERS=1;
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
			TRIMMING_MODE="$3";
			shift 3
		;;
		-rda|--run-de-novo)
			ASSEMBLY_FLAG=1;
			shift
		;;
		-rfa|--run-falcon)
			FALCON_FLAG=1;
			FALCON_MODE="$2";
			shift 2
		;;
		-enc|--encrypt)
			RUN_ENCRYPT=1;
			shift
		;;
		-dec|--decrypt)
			RUN_DECRYPT=1;
			shift
		;;
		-all|--run-all)
			TRIMMING_FLAG=1;
			TRIMMING_TYPE="TT";
			ASSEMBLY_FLAG=1;
			FALCON_FLAG=1;
			FALCON_MODE="SO";
			shift
		;;
		-*) # Unknown option
		echo -e "\033[1;34m[RFSC] \033[1;31m Invalid arg ($1)! \033[0m";
		echo -e "\033[1;34m[RFSC]\033[0m For more help, try: \033[0;33m./RFSC.sh -h \033[0m"
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

	echo "                                                                             "
	echo -e " \033[0;36m                       ____  _____ ____   ____                             \033[0m "
	echo -e " \033[0;36m                      |  _ \|  ___/ ___| / ___|                            \033[0m "
	echo -e " \033[0;36m                      | |_) | |_  \___ \| |                                \033[0m "
	echo -e " \033[0;36m                      |  _ <|  _|  ___) | |___                             \033[0m "
	echo -e " \033[0;36m                      |_| \_\_|   |____/ \____|                            \033[0m "
	echo "                                                                             "
	echo -e " \033[1;34m                           P I P E L I N E                                 \033[0m "
	echo "                                                                             "
	echo -e " \033[3;34m            A Reference-Free Sequence Classification Tool                  \033[0m "
	echo -e " \033[3;34m              for DNA sequences in metagenomic samples.                    \033[0m "
	echo "                                                                             "
	echo -e " \033[1;33m                      Usage: ./RFSC.sh [options]                           \033[0m "
	echo "                                                                             "
	echo "   -h,  --help            Show this help message and exit                    "
	echo "   -v,  --version         Show the version and some information              "
	echo "   -i,  --install         Installation of all the needed tools               "
	echo "                                                                             "
	echo -e "   -t,  --threads \033[0;34m<THREADS>\033[0m                                                  "
	echo "                          Number of threads to be used                       "
	echo "                                                                             "
	echo "   -bref, --build-ref-virus                                                  "
	echo "                          Build reference database of virus from NCBI        "
	echo "                                                                             "
	echo "   -gad,  --gen-adapters  Generate FASTA file with adapters                  "
	echo "                                                                             "
	echo -e "   -synt, --synthetic \033[0;34m[FILE1] : [FILE3]\033[0m                                      "
	echo "                          Generate a synthetical sequence using 3            "
	echo "                          reference files for testing purposes               "
	echo "                                                                             "
	echo -e "   -trim, --filter \033[0;34m<TOOL> <MODE>\033[0m                                             "
	echo "                          TOOL: Filter Reads using Trimmomatic (TT)          "
	echo "                                or using FASTP (FP)                          "
	echo "                          MODE: Paired End Reads (PE)                        "
	echo "                                or Single End Reads (SE)                     "
	echo "                                                                             "
	echo "   -rda, --run-de-novo    De-Novo Sequence Assembly                          "
	echo "                                                                             "
	echo -e "   -rfa, --run-falcon \033[0;34m<MODE>\033[0m                                                 "
	echo "                          Run Data Analysis with FALCON using only the       "
	echo "                          scaffolds (SO) or analysing by each                "
	echo "                          Read (RM)                                          "
	echo "                                                                             "
	echo "   -dec, --decrypt        Decrypt all files in /Data_Security/Decrypted_Data "
	echo "   -enc, --encrypt        Encrypt all files in /Data_Security/Encrypted_data "
	echo "                                                                             "
	echo "   -all, --run-all        Run all the options (considering real data)        "
	echo "                                                                             "
	exit 1;
fi
#
# ======================================================================
# VERSION
#
if [ "$SHOW_VERSION" -eq "1" ]; then

	echo "                                                         "
	echo "                          RFSC                           "
	echo "                                                         "
	echo "                      Version: 1.0                       "
	echo "                                                         "
	echo "                       IEETA/DETI                        "
	echo "             University of Aveiro, Portugal.             "
	echo "                                                         "
	exit 0;
fi	
#
# ======================================================================
# INSTALLATIONS
#
if [ "$INSTALL" -eq "1" ]; then
	./src/install_tools.sh
fi
#
# ======================================================================
# THREADS
#
if [ "$THREADS" -eq "1" ]; then
	THREADS_AVAILABLE=$GET_THREADS;
	echo -e "\033[1;34m[RFSC]\033[0m The system is now set to use $THREADS_AVAILABLE threads!"
fi
#
# ======================================================================
# BUILD REFERENCE VIRUS DATABASE
#
if [ "$BUILD_DB_VIRUS" -eq "1" ]; then
	./src/build_ref_db_virus.sh
fi	
#
# ======================================================================
# GENERATE FASTA ADAPTERS
#
if [ "$GEN_ADAPTERS" -eq "1" ]; then
	./src/gen_adapters.sh
fi	
#
# ===================================================================
#
if [[ "$GEN_SYNTHETIC" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Start Synthetic Sequence Generation!"
	GENERATE_SYNTHETIC "$REF_FILE1" "$REF_FILE2" "$REF_FILE3";
fi
#
# ===================================================================
#
if [[ "$TRIMMING_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Start Trimming the Sequence!"
	TRIMMING_SEQUENCE "$TRIMMING_TYPE";
fi
#
# ===================================================================
#
if [[ "$ASSEMBLY_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Start De-Novo Assembly!"
	SPADES_ASSEMBLY;
fi
#
# ===================================================================
#
if [[ "$FALCON_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Data Analysis!"
	FALCON_ANALYSIS "$FALCON_MODE";
fi
#
# ===================================================================
#
if [[ "$RUN_ENCRYPT" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Encryption!"
	ENCRYPT_DATA;
fi
#
# ===================================================================
#
if [[ "$RUN_DECRYPT" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Decryption!"
	DECRYPT_DATA;
fi
#