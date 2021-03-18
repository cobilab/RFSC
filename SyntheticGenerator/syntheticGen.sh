#!/bin/bash
#
# ============================================================== #
# =                                                            = #
# =                     Synthetic Generator                    = #
# =                                                            = #
# ============================================================== #
# =                                                            = #
# =                Usage: sh syntheticGen.sh [FILE]            = #
# =                                                            = #
# ============================================================== #
#
# ART ILLUMINA VALUES
RAND_SEED=0;				# -rs
SEQ_SYSTEM="HS25";			# -ss
ART_INPUT=$1;				# -i
ART_LEN=150;				# -l
ART_FOLD=40;				# -f
ART_SIZE=200;				# -m
ART_DSV=10;				# -s
ART_OUTPUT="SyntheticData/sample";	# -o

# GTO GENERATOR VALUES
GTO_NUM_SYMBOLS=5000;			# -n
GTO_LINE_SIZE=70;			# -l


if [ -n "$ART_INPUT" ]; then
	echo "Start FASTA Processing..."
else
	echo "Argument Error - Insert Input File"
	echo "Usage: bash syntheticGen.sh [FILE]"
	exit 0
fi

# CLEAR OLD GENERATED FILES
rm SyntheticData/sample*
#rm -rf SyntheticOutputs/out_spades_

# GENERATE RANDOM DNA GENOME
gto_genomic_gen_random_dna -n $GTO_NUM_SYMBOLS | gto_fasta_from_seq -n GTO_Genomic_Generator_Random_DNA_Sequence -l $GTO_LINE_SIZE > newSeq.fa

cat $ART_INPUT newSeq.fa > tmp.fa
echo "" >> tmp.fa
rm newSeq.fa $ART_INPUT
mv tmp.fa $ART_INPUT

# GENERATE SYNTHETIC SEQUENCE
art_illumina -rs $RAND_SEED -ss $SEQ_SYSTEM -sam -i $ART_INPUT -p -l $ART_LEN -f $ART_FOLD -m $ART_SIZE -s $ART_DSV -o $ART_OUTPUT

# GZIP FASTQ FILES
gzip SyntheticData/sample*.fq
