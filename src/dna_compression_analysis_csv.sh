#!/bin/bash

declare -a array_identifier
declare -a array_size
declare -a array_level
runs=0;
#files=0;
file_res=0;

printf "Identifier\tSize1\tSize2\tSize3\tSize4\tSize7\tSize12\tSize14\tLevel1\tLevel2\tLevel3\tLevel4\tLevel7\tLevel12\tLevel14\n" >> GeCo3_Output/Virus/geco3_Viral.csv

readarray -t array_identifier < GeCo3_Output/Virus/VIRAL_ID.txt

# for file in References/NCBI-Virus/*
# do
#     if [[ "$file" == "viral_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         array_identifier[$files]="${out_file%".fna.gz"}"
#         (( files++ ))
#     fi
# done

while read line; 
do
    if [[ "$runs" == "7" ]]; then
        printf "${array_identifier[$file_res]}\t${array_size[0]}\t${array_size[1]}\t${array_size[2]}\t${array_size[3]}\t${array_size[4]}\t${array_size[5]}\t${array_size[6]}\t${array_level[0]}\t${array_level[1]}\t${array_level[2]}\t${array_level[3]}\t${array_level[4]}\t${array_level[5]}\t${array_level[6]}\n" >> GeCo3_Output/Virus/geco3_Viral.csv
        (( file_res++ ))
        runs=0;
        array_size=()
        array_level=()
    fi

    array_size[$runs]=`echo $line | cut -d ' ' -f3`
    array_level[$runs]=`echo $line | cut -d ' ' -f8`
    array_level[$runs]=`echo "scale=5; ${array_level[$runs]}/2" | bc`
    (( runs++ ))
    
done <GeCo3_Output/Virus/geco_VIRAL.txt

# Input of the last set of data
printf "${array_identifier[$file_res]}\t${array_size[0]}\t${array_size[1]}\t${array_size[2]}\t${array_size[3]}\t${array_size[4]}\t${array_size[5]}\t${array_size[6]}\t${array_level[0]}\t${array_level[1]}\t${array_level[2]}\t${array_level[3]}\t${array_level[4]}\t${array_level[5]}\t${array_level[6]}\n" >> GeCo3_Output/Virus/geco3_Viral.csv

#echo "${array_size[0]}"
#echo "${array_level[6]}"