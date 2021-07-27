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
CLEAN=0;
#
BUILD_DB_VIRUS=0;
BUILD_DB_BACTERIA=0;
BUILD_DB_ARCHAEA=0;
BUILD_DB_PROTOZOA=0;
BUILD_DB_FUNGI=0;
BUILD_DB_PLANT=0;
#BUILD_DB_INVERTEBRATE=0;
#BUILD_DB_VERTEBRATE_MAMMALIAN=0;
#BUILD_DB_VERTEBRATE_OTHER=0;
BUILD_DB_MITOCHONDRIAL=0;
BUILD_DB_PLASTID=0;
#
GEN_ADAPTERS=0;
#
THREADS_AVAILABLE=`cat /proc/cpuinfo | awk '/^processor/{print $3}' | wc -l`;
THREADS=0;
GET_THREADS="";
#
SET_THRESHOLD_REF_BASED=0;
MAX_THRESHOLD_REF_BASED=70;
MIN_THRESHOLD_REF_BASED=1;
GET_MAX_THRESHOLD_REF_BASED="";
GET_MIN_THRESHOLD_REF_BASED="";
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
SET_LEN_COV=0;
NODE_LENGTH=100;
NODE_COVERAGE=3;
SET_NODE_LENGTH="";
SET_NODE_COVERAGE="";
#
FALCON_FLAG=0;
FALCON_MODE="";
BLASTN_REMOTE_FLAG=0;
#
RUN_DECRYPT=0;
RUN_ENCRYPT=0;
#
ORFFINDER_FLAG=0;
#
ORF_DATASETS=0;
ORF_TOOL="";
ORF_DOMAIN="";
#
NC_DNA_FLAG=0;
NC_DNA_DOMAIN="";
#
NC_AA_FLAG=0;
NC_AA_DOMAIN="";
#
GC_CONTENT_FLAG=0;
GC_CONTENT_DOMAIN="";
#
LEN_SEQ_FLAG=0;
LEN_SEQ_DOMAIN="";
#
TRAIN_TEST_DATASET_FLAG=0;
TRAIN_TEST_DATASET_PARTITION="";
#
LIMIT_SAMPLES_DATASET_FLAG=0;
LIMIT_SAMPLES_DATASET="";
#
RUN_GNB=0;
NUM_DOMAINS="";
PREDICTORS_CODE="";
#
RUN_KNN=0;
KNN_K="";
#
TEST_KNN_FLAG=0;
TEST_KNN_MODE="";
#
TEST_GNB_CV_FLAG=0;
TEST_GNB_DOMAIN="";
#
TEST_GNB_FLAG=0;
TEST_GNB_TRAIN_PERCENTAGE="";
#
ACCURACY_KNN_FLAG=0;
ACCURACY_KNN_MODE="";
ACCURACY_KNN_TRAIN="";
#
ACCURACY_GNB_FLAG=0;
ACCURACY_GNB_MODE="";
ACCURACY_GNB_TRAIN="";
ACCURACY_GNB_TRAIN_PERC="";
#
#
# ==================================================================
# VERIFICATION FUNCTIONS
#
CHECK_ADAPTERS() {
	if [ ! -f Input_Data/ReferenceBased/adapters.fa ]; then
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
	for file in Input_Data/ReferenceBased/*.fq.gz
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
		cp Input_Data/ReferenceBased/adapters.fa adapters.fa
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
# PARSE SCAFFOLDS INTO NODES (MINOR SCAFFOLDS)
#
PARSE_SCAFFOLDS() {
	echo -e "\033[1;34m[RFSC]\033[0m Start parsing scaffolds.fasta"
	mkdir GeneratedFiles/out_spades_/Nodes
	awk '/>/{filename="GeneratedFiles/out_spades_/Nodes/"NR".fasta"}; {print >filename}' GeneratedFiles/out_spades_/scaffolds.fasta
	echo -e "\033[1;34m[RFSC]\033[0m Parse finnished! The result can be find in GeneratedFiles/out_spades_/Nodes/"

	echo -e "\033[1;34m[RFSC]\033[0m Start filtering nodes based on Lenght & Coverage!"
	for file in GeneratedFiles/out_spades_/Nodes/*
	do
		REF=`grep ">" < $file`
		NODE_LEN=$(awk -F_ '{print $4}' <<< ${REF})
		NODE_COV=$(awk -F_ '{print $6}' <<< ${REF})
		COV_TO_INT=(${NODE_COV//./ })

		if [[ "$NODE_LEN" -gt "$NODE_LENGTH" ]] && [[ "${COV_TO_INT[0]}" -ge "$NODE_COVERAGE" ]]; then
			continue
		else
			rm $file
			echo -e "\033[1;34m[RFSC]\033[0m $file Node has been removed in the filtering process!"
		fi
	done
}
#
# ==================================================================
# FALCON ANALYSIS - SCAFFOLDS
#
FALCON_SO_MODE() {
	FALCON -n $THREADS_AVAILABLE -v -F -x Outputs/falcon_SO_results.txt GeneratedFiles/out_spades_/scaffolds.fasta References/NCBI-Virus/DB-viral.fa
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
FALCON_RM_MODE() {
	reads=0
	for file in GeneratedFiles/out_spades_/Nodes/*
	do
		array[ $reads ]=$file
		(( reads++ ))
	done
	#
	mkdir Outputs/FalconNodes
	len=${#array[@]}
	i=0
	while [[ $i != $len ]]
	do
		if [[ $(( ($len - $i) % 2)) == 0 ]]; then
			path1=${array[i]}
			path2=${array[i+1]}
			path_to_file1=(${path1//// })
			path_to_file2=(${path2//// })
			R1=${path_to_file1[3]}
			R2=${path_to_file2[3]}
			echo -e "\033[1;34m[RFSC]\033[0m Analysing Nodes $R1 & $R2 with $(($THREADS_AVAILABLE/2)) threads each!"
			FALCON -n $(($THREADS_AVAILABLE/2)) -v -F -x Outputs/FalconNodes/falcon_RM_"${R1}"_results.txt GeneratedFiles/out_spades_/Nodes/$R1 References/NCBI-Virus/DB-viral.fa | FALCON -n $(($THREADS_AVAILABLE/2)) -v -F -x Outputs/FalconNodes/falcon_RM_"${R2}"_results.txt GeneratedFiles/out_spades_/Nodes/$R2 References/NCBI-Virus/DB-viral.fa
			(( i+=2 ))
		else
			path=${array[i]}
			path_to_file=(${path//// })
			R=${path_to_file[3]}
			echo -e "\033[1;34m[RFSC]\033[0m Analysing Node $R with $(($THREADS_AVAILABLE)) threads!"
			FALCON -n $(($THREADS_AVAILABLE)) -v -F -x Outputs/FalconNodes/falcon_RM_"${R}"_results.txt GeneratedFiles/out_spades_/Nodes/$R References/NCBI-Virus/DB-viral.fa
			(( i++ ))
		fi 
	done
}
#
# ==================================================================
# SELECT RESULTS (AFTER FALCON)
#
SELECT_RESULTS() {
	echo -e "\033[1;34m[RFSC]\033[0m Starting selection procedure"
	mkdir Results/falcon_seq
	for file in Outputs/FalconNodes/*
	do
		readarray -t fasta_node <$file

		NLINES=${#fasta_node[@]}

		for (( i=0; i<$NLINES; i++ ));
		do
			PER=`echo "${fasta_node[i]}"|awk '{print $3}'`
			FIND_MATCH=`echo "$PER > $MAX_THRESHOLD_REF_BASED" | bc`
			SECOND_PHASE=`echo "$PER > $MIN_THRESHOLD_REF_BASED && $PER < $MAX_THRESHOLD_REF_BASED" | bc`

			file=${file#"Outputs/FalconNodes/falcon_RM_"}
			file=${file%"_results.txt"}

			if [[ $FIND_MATCH -eq "1" ]]; then
				GENOME=`echo "${fasta_node[i]}"|awk '{print $4}'`

				echo -e "\033[1;34m[RFSC]\033[0m Moving $file to Results/falcon_seq"
				mv GeneratedFiles/out_spades_/Nodes/$file Results/falcon_seq

				ALREADY_FOUND=`grep -x $GENOME Results/ref_result.txt`
				if [[ $ALREADY_FOUND != "" ]]; then
					echo -e "\033[1;34m[RFSC]\033[0m Sequence already stated in Results/ref_result.txt"
				else
					printf "$GENOME\n" >> Results/ref_result.txt
				fi

				break
			elif [[ $SECOND_PHASE -eq "1" ]]; then
				echo -e "\033[1;34m[RFSC]\033[0m Moving $file to Input_Data/ReferenceFree"
				mv GeneratedFiles/out_spades_/Nodes/$file Input_Data/ReferenceFree
				break
			else
				echo -e "\033[1;34m[RFSC]\033[0m Deleting $file"
				rm GeneratedFiles/out_spades_/Nodes/$file
				break
			fi

		done
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
		PARSE_SCAFFOLDS;
		FALCON_RM_MODE;
		SELECT_RESULTS;
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
# BLASTN ANALYSIS
#
BLASTN_ANALYSIS() {
	PARSE_SCAFFOLDS;
	mkdir Outputs/BlastnNodes
	for file in GeneratedFiles/out_spades_/Nodes/*
	do
		file_name=$(basename $file);
		echo -e "\033[1;34m[RFSC]\033[0m Blastn is now processing $file"
		blastn -db nt -task blastn-short -query $file -remote > Outputs/BlastnNodes/$file_name.txt
	done
}
#
# ==================================================================
# ORF FINDER SEARCH (ORFfinder)
#
ORF_SEARCH() {
	mkdir ORFs/NodesProteins
	for file in Input_Data/ReferenceFree/*
	do
		node_file=$(basename $file);
		echo -e "\033[1;34m[RFSC]\033[0m ORFfinder is now processing $node_file"
		node_file=${node_file%".fasta"}
		./ORFs/ORFfinder -s 1 -in Input_Data/ReferenceFree/$node_file.fasta -outfmt 0 -out ORFs/NodesProteins/$node_file.protein.fasta
	done
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
		-tmm|--set-threshold-max-min)
			SET_THRESHOLD_REF_BASED=1;
			GET_MAX_THRESHOLD_REF_BASED="$2";
			GET_MIN_THRESHOLD_REF_BASED="$3";
			shift 3
		;;
		-dlc|--set-len-cov)
			SET_LEN_COV=1;
			SET_NODE_LENGTH="$2";
			SET_NODE_COVERAGE="$3";
			shift 3
		;;
		-bviral|--build-ref-virus)
			BUILD_DB_VIRUS=1;
			shift
		;;
		-bbact|--build-ref-bacteria)
			BUILD_DB_BACTERIA=1;
			shift
		;;
		-barch|--build-ref-archaea)
			BUILD_DB_ARCHAEA=1;
			shift
		;;
		-bprot|--build-ref-protozoa)
			BUILD_DB_PROTOZOA=1;
			shift
		;;
		-bfung|--build-ref-fungi)
			BUILD_DB_FUNGI=1;
			shift
		;;
		-bplan|--build-ref-plant)
			BUILD_DB_PLANT=1;
			shift
		;;
		-bmito|--build-ref-mitochondrial)
			BUILD_DB_MITOCHONDRIAL=1;
			shift
		;;
		-bplas|--build-ref-plastid)
			BUILD_DB_PLASTID=1;
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
		-rbr|--run-blastn-remote)
			BLASTN_REMOTE_FLAG=1;
			shift
		;;
		-orf|--orf-finder)
			ORFFINDER_FLAG=1;
			shift
		;;
		-orfd|--orf-dataset)
			ORF_DATASETS=1;
			ORF_TOOL="$2";
			ORF_DOMAIN="$3";
			shift 3
		;;
		-enc|--encrypt)
			RUN_ENCRYPT=1;
			shift
		;;
		-dec|--decrypt)
			RUN_DECRYPT=1;
			shift
		;;
		-ball|--build-ref-all)
			BUILD_DB_VIRUS=1;
			BUILD_DB_BACTERIA=1;
			BUILD_DB_ARCHAEA=1;
			BUILD_DB_PROTOZOA=1;
			BUILD_DB_FUNGI=1;
			BUILD_DB_PLANT=1;
			BUILD_DB_INVERTEBRATE=1;
			BUILD_DB_VERTEBRATE_MAMMALIAN=1;
			BUILD_DB_VERTEBRATE_OTHER=1;
			BUILD_DB_MITOCHONDRIAL=1;
			BUILD_DB_PLASTID=1;
			shift
		;;
		-ncd|--nc-dna-csv)
			NC_DNA_FLAG=1;
			NC_DNA_DOMAIN="$2";
			shift 2
		;;
		-nca|--nc-aa-csv)
			NC_AA_FLAG=1;
			NC_AA_DOMAIN="$2";
			shift 2
		;;
		-gc|--gc-content-csv)
			GC_CONTENT_FLAG=1;
			GC_CONTENT_DOMAIN="$2";
			shift 2
		;;
		-lenseq|--len-dna-aa-csv)
			LEN_SEQ_FLAG=1;
			LEN_SEQ_DOMAIN="$2";
			shift 2
		;;
		-train-test|--train-test-dataset-csv)
			TRAIN_TEST_DATASET_FLAG=1;
			TRAIN_TEST_DATASET_PARTITION="$2";
			shift 2
		;;
		-sdataset|--small-dataset-csv)
			LIMIT_SAMPLES_DATASET_FLAG=1;
			LIMIT_SAMPLES_DATASET="$2";
			shift 2
		;;
		-gnb|--run-gaussian-naive-bayes-classifier)
			RUN_GNB=1;
			NUM_DOMAINS="$2";
			PREDICTORS_CODE="$3";
			shift 3
		;;
		-knn|--run-k-nearest-neighbor-classifier)
			RUN_KNN=1;
			KNN_K="$2";
			shift 2
		;;
		-testKNN)
			TEST_KNN_FLAG=1;
			TEST_KNN_MODE="$2";
			shift 2
		;;
		-testGNB)
			TEST_GNB_FLAG=1;
			TEST_GNB_TRAIN_PERCENTAGE="$2";
			shift 2
		;;
		-testGNB-CV|--testGNB-CrossV)
			TEST_GNB_CV_FLAG=1;
			TEST_GNB_DOMAIN="$2";
			shift 2
		;;
		-aKNN|--accuracy-KNN)
			ACCURACY_KNN_FLAG=1;
			ACCURACY_KNN_MODE="$2";
			ACCURACY_KNN_TRAIN="$3";
			shift 3
		;;
		-aGNB|--accuracy-GNB)
			ACCURACY_GNB_FLAG=1;
			ACCURACY_GNB_MODE="$2";
			ACCURACY_GNB_TRAIN="$3";
			ACCURACY_GNB_TRAIN_PERC="$4"
			shift 3
		;;
		-clc|--clean-all)
			CLEAN=1;
			shift
		;;
		-all|--run-all)
			CLEAN=1;
			SET_THRESHOLD_REF_BASED=1;
			GET_MAX_THRESHOLD_REF_BASED="70";
			GET_MIN_THRESHOLD_REF_BASED="1";
			SET_LEN_COV=1;
			SET_NODE_LENGTH="100";
			SET_NODE_COVERAGE="3";
			TRIMMING_FLAG=1;
			TRIMMING_TYPE="TT";
			TRIMMING_MODE="PE";
			ASSEMBLY_FLAG=1;
			FALCON_FLAG=1;
			FALCON_MODE="RM";
			RUN_KNN=1;
			KNN_K="2";
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
	echo "   -dec, --decrypt        Decrypt all files in /Data_Security/Decrypted_Data "
	echo "   -enc, --encrypt        Encrypt all files in /Data_Security/Encrypted_data "
	echo "                                                                             "
	echo -e "   -tmm,  --set-threshold-max-min \033[0;34m<MAX> <MIN>\033[0m                                "
	echo "                          Set Max & Min thresholds for percentage            "
	echo "                          similarity in reference based analysis             "
	echo "                                                                             "
	echo -e "   -dlc,  --set-len-cov \033[0;34m<LEN> <COV>\033[0m                                          "
	echo "                          Define the Length and Coverage values              "
	echo "                          for the scaffolds filtering process                "
	echo "                                                                             "
	echo "   -gad,  --gen-adapters  Generate FASTA file with adapters                  "
	echo "                                                                             "
	echo -e "   -synt, --synthetic \033[0;34m[FILE1]:[FILE3]\033[0m                                        "
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
	echo "   -orf, --orf-finder     Perform DNA sequence translation for amino acids,  "
	echo "                          finds all open reading frames (ORF) and remove     "
	echo "                          stop codons                                        "
	echo "                                                                             "
	echo -e "   -orfd, --orf-dataset \033[0;34m<TOOL> <DOMAIN>                                     \033[0m "
	echo "                          Converts nucleotide sequences presented in the     "
	echo "                          NCBI databasesinto protein sequences               "
	echo "                          TOOL: orfM (-orfm) [DEFAULT]                       "
	echo "                                or using ORF-Finder (-orffinder)             "
	echo "                          DOMAIN: --viral, --bacteria, --archaea, ...        "
	echo "                                                                             "
	echo -e "   -ncd, --nc-dna-csv \033[0;34m<DOMAIN>                                              \033[0m "
	echo "                          Compresses and generates a CSV file for the        "
	echo "                          DNA NCBI datasets                                  "
	echo "                          DOMAIN: --viral, --bacteria, --archaea, ...        "
	echo "                                                                             "
	echo -e "   -nca, --nc-aa-csv \033[0;34m<DOMAIN>                                               \033[0m "
	echo "                          Compresses and generates a CSV file for the        "
	echo "                          AA NCBI datasets                                   "
	echo "                          DOMAIN: --viral, --bacteria, --archaea, ...        "
	echo "                                                                             "
	echo -e "   -gc, --gc-content-csv \033[0;34m<DOMAIN>                                           \033[0m "
	echo "                          Analyses the percentage of GC-Content in each      "
	echo "                          sequece of the choosen NCBI database               "
	echo "                          DOMAIN: --viral, --bacteria, --archaea, ...        "
	echo "                                                                             "
	echo -e "   -lenseq, --len-dna-aa-csv \033[0;34m<DOMAIN>                                       \033[0m "
	echo "                          Analyses the Length (DNA & AA) sequences the       "
	echo "                          choosen NCBI database                              "
	echo "                          DOMAIN: --viral, --bacteria, --archaea, ...        "
	echo "                                                                             "
	echo -e " \033[1;33m                       B U I L D    D A T A B A S E S                      \033[0m "
	echo "                                                                             "
	echo "   -bviral, --build-ref-virus                                                "
	echo -e "                          Build reference database for \033[1;36mvirus\033[0m from NCBI       "
	echo "   -bbact,  --build-ref-bacteria                                             "
	echo -e "                          Build reference database for \033[1;36mbacterias\033[0m from NCBI   "
	echo "   -barch,  --build-ref-archaea                                              "
	echo -e "                          Build reference database for \033[1;36marchaeas\033[0m from NCBI    "
	echo "   -bprot,  --build-ref-protozoa                                             "
	echo -e "                          Build reference database for \033[1;36mprotozoa\033[0m from NCBI    "
	echo "   -bfung,  --build-ref-fungi                                                "
	echo -e "                          Build reference database for \033[1;36mfungi\033[0m from NCBI       "
	echo "   -bplan,  --build-ref-plant                                                "
	echo -e "                          Build reference database for \033[1;36mplant\033[0m from NCBI       "
	echo "   -bmito,  --build-ref-mitochondrial                                        "
	echo -e "                          Build reference database for \033[1;36mmitochondrial\033[0m from NCBI"
	echo "   -bplas,  --build-ref-plastid                                              "
	echo -e "                          Build reference database for \033[1;36mplastid\033[0m from NCBI     "
	echo -e " \033[1;33m                - - - - - - - - - - - - - - - - - - - - - -                \033[0m "
	echo "                                                                             "
	echo -e " \033[1;33m              R E F E R E N C E   B A S E D   A P P R O A C H              \033[0m "
	echo "                                                                             "
	echo -e "   -rfa, --run-falcon \033[0;34m<MODE>\033[0m                                                 "
	echo "                          Run Data Analysis with FALCON using only the       "
	echo "                          scaffolds (SO) or analysing by each                "
	echo "                          Read (RM)                                          "
	echo "                                                                             "
	echo "   -rbr, --run-blastn-remote                                                 "
	echo "                          Run Data Analysus with Blast+ using remote         "
	echo "                          access to NCBI databases                           "
	echo -e " \033[1;33m                - - - - - - - - - - - - - - - - - - - - - -                \033[0m "
	echo "                                                                             "
	echo -e "   -train-test, --train-test-dataset-csv \033[0;34m<TRAIN_PARTITION>                  \033[0m "
	echo "                          Divide the dataset into a train and test dataset   "
	echo "                          (Usefull for testing the KNN Classifier)           "
	echo "                          TRAIN_PARTITION: value (0..1) for the train partition"
	echo "                                                                             "
	echo -e "   -sdataset, --small-dataset-csv \033[0;34m<MAX_SAMPLES>                             \033[0m "
	echo "                          Create a small dataset with a maximum of samples   "
	echo "                          for each domain                                    "
	echo "                          MAX_SAMPLES: Maximum samples per domain            "
	echo "                                                                             "
	echo -e " \033[1;33m               R E F E R E N C E   F R E E   A P P R O A C H               \033[0m "
	echo "                                                                             "
	echo -e "   -gnb, --run-gaussian-naive-bayes-classifier \033[0;34m<NUM_DOMAINS> <PREDICTORS>\033[0m    "
	echo "                          NUM_DOMAINS: Number of domains supported           "
	echo "                          PREDICTORS: Code regarding the desired predictors: "
	echo "                '1111' -> All predictors                                     "
	echo "                '0001' -> Only Nucleotide Compression                        "
	echo "                '0010' -> Only AA Compression                                "
	echo "                '0011' -> Only GC-Content                                    "
	echo "                '0100' -> Only DNA Length                                    "
	echo "                '0101' -> Only AA Length                                     "
	echo "                '0110' -> DNA & AA Compression                               "
	echo "                '0111' -> DNA & AA Compression + GC-Content                  "
	echo "                '1000' -> DNA & AA Compression + GC-Content + DNA Length     "
	echo "                '1001' -> DNA & AA Compression + GC-Content + DNA & AA Length"
	echo "                '1010' -> DNA & AA Compression + DNA & AA Length             "
	echo "                                                                             "
	echo -e "   -knn, --run-k-nearest-neighbor-classifier \033[0;34m<K>                            \033[0m "
	echo "                          K: Number of neighbors for the classifier          "
	echo -e " \033[1;33m                - - - - - - - - - - - - - - - - - - - - - -                \033[0m "
	echo "                                                                             "
	echo -e " \033[1;33m         T E S T    C L A S S I F I E R S    P E R F O R M A N C E    \033[0m      "
	echo "                                                                             "
	echo -e "   -testKNN \033[0;34m<MODE>                                                          \033[0m "
	echo "                          Deep test the KNN Classifier                       "
	echo "                          MODE: Test Database against Train Database (--test)"
	echo "                                Cross-Validation (--viral, --bacteria, ...)  "
	echo "                                                                             "
	echo -e "   -testGNB \033[0;34m<PERCENTAGE>                                                    \033[0m "
	echo "                          Deep test the GNB Classifier with Test & Train Datasets "
	echo "                          PERCENTAGE: Value between 0..1 representing the    "
	echo "                                      percentage of the dataset reserved for "
	echo "                                      the training                           "
	echo "                                                                             "
	echo -e "   -testGNB-CV, --testGNB-CrossV \033[0;34m<DOMAIN>                                   \033[0m "
	echo "                          Deep test the GNB Classifier with Cross-Validation "
	echo "                          DOMAIN: --viral, --bacteria, --archaea, ...        "
	echo "                                                                             "
	echo -e "   -aKNN, --accuracy-KNN \033[0;34m<AC-MODE> <T-MODE>                                 \033[0m "
	echo "                          Analyse the accuracy of the KNN Classifier when    "
	echo "                          testing it against a known dataset                 "
	echo "                          AC-MODE: Simple Accuracy Mean (Accuracy)           "
	echo "                                   Weighted F1-Score (F1Score)               "
	echo "                          T-MODE: Cross-Validation Method (CV)               "
	echo "                                  Train-Test Database (Test)                 "
	echo "                                                                             "
	echo -e "   -aGNB, --accuracy-GNB \033[0;34m<AC-MODE> <T-MODE> <TRAIN-PERCENTAGE>              \033[0m "
	echo "                          Analyse the accuracy of the KNN Classifier when    "
	echo "                          testing it against a known dataset                 "
	echo "                          AC-MODE: Simple Accuracy Mean (Accuracy)           "
	echo "                                   Weighted F1-Score (F1Score)               "
	echo "                          T-MODE: Cross-Validation Method (CV)               "
	echo "                                  Train-Test Database (Test)                 "
	echo "                          TRAIN-PERCENTAGE: Percentage of the dataset used   "
	echo "                                            for training purpouses           "
	echo -e " \033[1;33m                - - - - - - - - - - - - - - - - - - - - - -                \033[0m "
	echo "                                                                             "
	echo "   -clc, --clean-all      Clean all generated files (Including Results)      "
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
# REFERENCE BASED PERCENTAGE (MAX MIN) THRESHOLDS
#
if [ "$SET_THRESHOLD_REF_BASED" -eq "1" ]; then
	MAX_THRESHOLD_REF_BASED=$GET_MAX_THRESHOLD_REF_BASED;
	MIN_THRESHOLD_REF_BASED=$GET_MIN_THRESHOLD_REF_BASED;
	echo -e "\033[1;34m[RFSC]\033[0m The system is now set to use $MAX_THRESHOLD_REF_BASED & $MIN_THRESHOLD_REF_BASED as thresholds for reference based analysis!"
fi
#
# ======================================================================
# SET LENGTH & COVERAGE VALUES FOR SCAFFOLDS FILTERING
#
if [ "$SET_LEN_COV" -eq "1" ]; then
	NODE_LENGTH=$SET_NODE_LENGTH;
	NODE_COVERAGE=$SET_NODE_COVERAGE;
	echo -e "\033[1;34m[RFSC]\033[0m The Coverage value for each node is now set to $SET_NODE_COVERAGE & length to $SET_NODE_LENGTH!"
fi
#
# ======================================================================
# BUILD REFERENCE VIRAL DATABASE
#
if [ "$BUILD_DB_VIRUS" -eq "1" ]; then
	cd References/NCBI-Virus/
	echo -e "\033[1;34m[RFSC]\033[0m Building viral database at References/NCBI-Virus/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --viral
	gunzip DB-viral.fa.gz
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE BACTERIAL DATABASE
#
if [ "$BUILD_DB_BACTERIA" -eq "1" ]; then
	cd References/NCBI-Bacteria/
	echo -e "\033[1;34m[RFSC]\033[0m Building bacterias database at References/NCBI-Archaea/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --bacteria
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE ARCHAEAS DATABASE
#
if [ "$BUILD_DB_ARCHAEA" -eq "1" ]; then
	cd References/NCBI-Archaea/
	echo -e "\033[1;34m[RFSC]\033[0m Building archaeas database at References/NCBI-Archaea/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --archaea
	gunzip DB-archaea.fa.gz
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE PROTOZOA DATABASE
#
if [ "$BUILD_DB_PROTOZOA" -eq "1" ]; then
	cd References/NCBI-Protozoa/
	echo -e "\033[1;34m[RFSC]\033[0m Building protozoa database at References/NCBI-Protozoa/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --protozoa
	gunzip DB-protozoa.fa.gz
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE FUNGI DATABASE
#
if [ "$BUILD_DB_FUNGI" -eq "1" ]; then
	cd References/NCBI-Fungi/
	echo -e "\033[1;34m[RFSC]\033[0m Building fungi database at References/NCBI-Fungi/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --fungi
	gunzip DB-fungi.fa.gz
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE PLANT DATABASE
#
if [ "$BUILD_DB_PLANT" -eq "1" ]; then
	cd References/NCBI-Plant/
	echo -e "\033[1;34m[RFSC]\033[0m Building plant database at References/NCBI-Plant/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --plant
	gunzip DB-plant.fa.gz
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE MITOCHONDRIAL DATABASE
#
if [ "$BUILD_DB_MITOCHONDRIAL" -eq "1" ]; then
	cd References/NCBI-Mitochondrial/
	echo -e "\033[1;34m[RFSC]\033[0m Building mitochondrial database at References/NCBI-Mitochondrial/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --mitochondrion
	gunzip DB-mitochondrion.fa.gz
	cd ../..
fi
#
# ======================================================================
# BUILD REFERENCE PLASTID DATABASE
#
if [ "$BUILD_DB_PLASTID" -eq "1" ]; then
	cd References/NCBI-Plastid/
	echo -e "\033[1;34m[RFSC]\033[0m Building plastid database at References/NCBI-Plastid/";
	./../../src/BUILD_DB.sh --threads $THREADS_AVAILABLE --plastid
	gunzip DB-plastid.fa.gz
	cd ../..
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
	echo -e "\033[1;34m[RFSC]\033[0m Starting Data Analysis with FALCON!"
	FALCON_ANALYSIS "$FALCON_MODE";
fi
#
# ===================================================================
#
if [[ "$BLASTN_REMOTE_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Data Analysis with Remote Blastn!"
	BLASTN_ANALYSIS;
fi
#
# ===================================================================
#
if [[ "$ORFFINDER_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Open Reading Frames (ORFs) Search!"
	ORF_SEARCH;
fi
#
# ===================================================================
#
if [[ "$ORF_DATASETS" -eq "1" ]]; then
	./src/protein_conversion_DB.sh $ORF_TOOL $ORF_DOMAIN
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
if [[ "$CLEAN" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Cleaning all files..."
	read -p "Are you sure you want to proceed? " -n 1 -r
	echo    # (optional) move to a new line
	if [[ $REPLY =~ ^[Yy]$ ]]; then
    	./clean.sh
	fi
fi
#
# ===================================================================
# ==== R E F E R E N C E   F R E E   C L A S S I F I C A T I O N ====
# ===================================================================
#
if [[ "$NC_DNA_FLAG" -eq "1" ]]; then
	./src/dna_compression.sh $NC_DNA_DOMAIN
fi
#
if [[ "$NC_AA_FLAG" -eq "1" ]]; then
	./src/protein_compression.sh $NC_AA_DOMAIN
fi
#
if [[ "$GC_CONTENT_FLAG" -eq "1" ]]; then
	./src/GCcontent_DB_analysis.sh $GC_CONTENT_DOMAIN
fi
#
if [[ "$LEN_SEQ_FLAG" -eq "1" ]]; then
	./src/length_DB_analysis.sh $LEN_SEQ_DOMAIN
fi
#
if [[ "$TRAIN_TEST_DATASET_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Generating the Train and Test Dataset"
	./src/train_test_dataset_generator.sh $TRAIN_TEST_DATASET_PARTITION
fi
#
if [[ "$LIMIT_SAMPLES_DATASET_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Generating the Dataset"
	./src/all_in_one_csv_generator.sh $LIMIT_SAMPLES_DATASET
fi
#
if [[ "$RUN_GNB" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Processing Input Sequences"
	./src/get_input_predictors.sh
	./src/ref_free_analysis.sh GNB 0 $NUM_DOMAINS $PREDICTORS_CODE
fi
#
if [[ "$RUN_KNN" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Processing Input Sequences!"
	./src/get_input_predictors.sh
	./src/ref_free_analysis.sh KNN $KNN_K 0 NULL
fi
#
if [[ "$TEST_KNN_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Testing the KNN Classifier!"
	./src/deepTest_KNN.sh $TEST_KNN_MODE
fi
#
if [[ "$TEST_GNB_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Testing the GNB Classifier using $TEST_GNB_TRAIN_PERCENTAGE % of the dataset for training!"
	./src/deepTest_GNB_TT.sh $TEST_GNB_TRAIN_PERCENTAGE --viral --bacteria --archaea --fungi --plant --protozoa --mitochondrion --plastid
fi
#
if [[ "$TEST_GNB_CV_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Testing the GNB Classifier with 5-Fold Cross-Validation!"
	./src/deepTest_GNB.sh $TEST_GNB_DOMAIN
fi
#
if [[ "$ACCURACY_KNN_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Analysing the KNN Classifier Accuracy using $ACCURACY_KNN_MODE!"
	if [[ $ACCURACY_KNN_TRAIN == "Test" ]]; then
		python3 Tests/accuracy_KNN.py $ACCURACY_KNN_MODE
	elif [[ $ACCURACY_KNN_TRAIN == "CV" ]]; then
		python3 Tests/accuracy_CV.py KNN/CrossValidation $ACCURACY_KNN_MODE
	else
		echo -e "\033[1;34m[RFSC]\033[0m $ACCURACY_KNN_MODE Testing mode is not recognized. Please insert a valid one!"
	fi
fi
#
if [[ "$ACCURACY_GNB_FLAG" -eq "1" ]]; then
	echo -e "\033[1;34m[RFSC]\033[0m Starting Analysing the GNB Classifier Accuracy using $ACCURACY_GNB_MODE!"
	if [[ $ACCURACY_GNB_TRAIN == "Test" ]]; then
		python3 Tests/accuracy_GNB.py GNB/CrossValidation $ACCURACY_GNB_TRAIN_PERC $ACCURACY_GNB_MODE
	elif [[ $ACCURACY_GNB_TRAIN == "CV" ]]; then
		python3 Tests/accuracy_CV.py GNB/CrossValidation $ACCURACY_GNB_MODE
	else
		echo -e "\033[1;34m[RFSC]\033[0m $ACCURACY_KNN_MODE Testing mode is not recognized. Please insert a valid one!"
	fi
fi
#