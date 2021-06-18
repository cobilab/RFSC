#!/bin/bash
#
# run: ./src/length_sequence_csv.sh Virus VIRAL length_Viral.csv
#      ./src/length_sequence_csv.sh Bacteria BACTERIA length_Bacteria.csv
#      ./src/length_sequence_csv.sh Archaea ARCHAEA length_Archaea.csv
#      ./src/length_sequence_csv.sh Protozoa PROTOZOA length_Protozoa.csv
#      ./src/length_sequence_csv.sh Fungi FUNGI length_Fungi.csv
#      ./src/length_sequence_csv.sh Plant PLANT length_Plant.csv
#      ./src/length_sequence_csv.sh Mitochondrial MITOCHONDRIAL length_Mitochondrial.csv
#      ./src/length_sequence_csv.sh Plastid PLASTID length_Plastid.csv
#
declare -a array_identifier
declare -a array_dna_symbols
declare -a array_aa_symbols

file_res=0;

printf "Identifier\tDNA_Length\tAA_Length\n" >> Analysis/LengthSeq/$1/$3

readarray -t array_identifier < Analysis/LengthSeq/$1/$2_DNA_ID.txt
readarray -t array_dna_symbols < Analysis/LengthSeq/$1/Length_DNA_$2.txt
readarray -t array_aa_symbols < Analysis/LengthSeq/$1/Length_AA_$2.txt

while read line; 
do

    printf "${array_identifier[$file_res]}\t${array_dna_symbols[$file_res]}\t${array_aa_symbols[$file_res]}\n" >> Analysis/LengthSeq/$1/$3
    (( file_res++ ))

done <Analysis/LengthSeq/$1/$2_DNA_ID.txt