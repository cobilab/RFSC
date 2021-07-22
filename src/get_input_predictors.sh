#!/bin/bash

#
MAX_LEN_DNA_DB=`python3 src/find_greatest_length.py DNA`
MAX_LEN_AA_DB=`python3 src/find_greatest_length.py AA`
#
if [ -f Input_Data/ReferenceFree/*.gz ]; then
    echo "Start processing Sequence(s)"
elif [ -f Input_Data/ReferenceFree/*.fna ]; then
    for file in Input_Data/ReferenceFree/*.fna
    do
        gzip $file
    done
elif [ -f Input_Data/ReferenceFree/*.fasta ]; then
    for file in Input_Data/ReferenceFree/*.fasta
    do
        gzip $file
    done
else
    echo -e "\033[1;34m[RFSC]\033[0m No .gz, .fna or .fasta files found in Input_Data/ReferenceFree to process. \n\033[1;34m[RFSC] !!\033[0m Please insert your sequence(s) in this directory!"
    exit 0
fi

# Create dir for Generated ORFs (From input sequences)
if [ ! -d Input_Data/ReferenceFree/GeneratedORFs ]; then
    mkdir Input_Data/ReferenceFree/GeneratedORFs
fi

# Create dir for Generated Predictor Values (From input sequences)
if [ ! -d Input_Data/ReferenceFree/GeneratedPredictorValues ]; then
    mkdir Input_Data/ReferenceFree/GeneratedPredictorValues
fi

# Create CSV for Input Sequences
printf "Identifier,GeCo3_1,AC_1,GeCo3_3,AC_3,GeCo3_7,AC_7,GC,DNA_LEN,AA_LEN\n" >> Input_Data/ReferenceFree/GeneratedPredictorValues/Input_Predictors.csv

for file in Input_Data/ReferenceFree/*.gz
do
    base_file=$(basename $file);
    base_file="${base_file%".gz"}"

    # ORF
    echo -e "\033[1;34m[RFSC]\033[0m ORFm is now processing $base_file"
    ./ORFs/orfm/orfm $file > Input_Data/ReferenceFree/GeneratedORFs/PROTEIN.fasta
    cat Input_Data/ReferenceFree/GeneratedORFs/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > Input_Data/ReferenceFree/GeneratedORFs/$base_file.fna
    rm Input_Data/ReferenceFree/GeneratedORFs/PROTEIN.fasta

    echo -e "\033[1;34m[RFSC]\033[0m Getting predictors for $base_file: "

    # GC-Content
    echo -e "\033[1;34m[RFSC]\033[0m Analysing GC-Content $base_file"
    GC=`./src/GC_analysis.sh $file ZIPPED`
    GC_CONTENT=`echo $GC | cut -d ' ' -f6`

    # DNA Compression
    echo -e "\033[1;34m[RFSC]\033[0m Analysing GeCo3 Compression $base_file"

    zcat $file | grep -v ">" | tr -d -c "ACGT" > Input_Data/ReferenceFree/GeneratedPredictorValues/TO_COMPRESS
    GECO_OUTPUT1=`GeCo3 -v -l 1 Input_Data/ReferenceFree/GeneratedPredictorValues/TO_COMPRESS | sed '1,6d' | sed '2d'`
    GECO_OUTPUT3=`GeCo3 -v -l 3 Input_Data/ReferenceFree/GeneratedPredictorValues/TO_COMPRESS | sed '1,6d' | sed '2d'`
    GECO_OUTPUT7=`GeCo3 -v -l 7 Input_Data/ReferenceFree/GeneratedPredictorValues/TO_COMPRESS | sed '1,6d' | sed '2d'`

    GECO_VALUE1=`echo $GECO_OUTPUT1 | cut -d ' ' -f8`
    GECO_VALUE3=`echo $GECO_OUTPUT3 | cut -d ' ' -f8`
    GECO_VALUE7=`echo $GECO_OUTPUT7 | cut -d ' ' -f8`

    GECO3_1=`echo "scale=5; ${GECO_VALUE1}/2" | bc`
    GECO3_3=`echo "scale=5; ${GECO_VALUE3}/2" | bc`
    GECO3_7=`echo "scale=5; ${GECO_VALUE7}/2" | bc`
    rm Input_Data/ReferenceFree/GeneratedPredictorValues/TO_COMPRESS*

    # AA Compression
    echo -e "\033[1;34m[RFSC]\033[0m Analysing AC Compression $base_file"
    AC_OUTPUT1=`AC -v -l 1 Input_Data/ReferenceFree/GeneratedORFs/$base_file.fna | sed '3,4d' | sed 'N;s/\n/ /'`
    AC_OUTPUT3=`AC -v -l 3 Input_Data/ReferenceFree/GeneratedORFs/$base_file.fna | sed '3,4d' | sed 'N;s/\n/ /'`
    AC_OUTPUT7=`AC -v -l 7 Input_Data/ReferenceFree/GeneratedORFs/$base_file.fna | sed '3,4d' | sed 'N;s/\n/ /'`

    AC_1=`echo $AC_OUTPUT1 | cut -d ' ' -f16`
    AC_3=`echo $AC_OUTPUT3 | cut -d ' ' -f16`
    AC_7=`echo $AC_OUTPUT7 | cut -d ' ' -f16`

    #echo $AC_OUTPUT | cut -d ' ' -f16 >> Input_Data/ReferenceFree/GeneratedPredictorValues/AC.txt
    rm Input_Data/ReferenceFree/GeneratedORFs/$base_file.fna.co

    # DNA Length Analysis
    echo -e "\033[1;34m[RFSC]\033[0m Analysing Nucleotide Length $base_file"
    DNA_LENGTH=`./src/length_sequence.sh $file DNA ZIPPED`
    DNA_LENGTH=`echo "scale=16; ${DNA_LENGTH}/${MAX_LEN_DNA_DB}" | bc`

    # AA Length Analysis
    echo -e "\033[1;34m[RFSC]\033[0m Analysing Amino-Acid Length $base_file"
    AA_LENGTH=`./src/length_sequence.sh Input_Data/ReferenceFree/GeneratedORFs/$base_file.fna AA UNZIPPED`
    AA_LENGTH=`echo "scale=16; ${AA_LENGTH}/${MAX_LEN_AA_DB}" | bc`

    printf "$base_file,0$GECO3_1,$AC_1,0$GECO3_3,$AC_3,0$GECO3_7,$AC_7,0$GC_CONTENT,0$DNA_LENGTH,0$AA_LENGTH\n" >> Input_Data/ReferenceFree/GeneratedPredictorValues/Input_Predictors.csv
done

rm -rf Input_Data/ReferenceFree/GeneratedORFs