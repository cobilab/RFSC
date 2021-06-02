#!/bin/bash

# run: ./src/compression_analysis_csv.sh Virus VIRAL geco3_Viral.csv NM
#      ./src/compression_analysis_csv.sh Fungi FUNGI geco3_Fungi.csv NM
#      ./src/compression_analysis_csv.sh Archaea ARCHAEA geco3_Archaea.csv NM
#      ./src/compression_analysis_csv.sh Plant PLANT geco3_Plant.csv NM
#      ./src/compression_analysis_csv.sh Protozoa PROTOZOA geco3_Protozoa.csv NM
#      ./src/compression_analysis_csv.sh Plastid PLASTID geco3_Plastid.csv NM
#      ./src/compression_analysis_csv.sh Mitochondrial MITOCHONDRIAL geco3_Mitochondrial.csv NM


declare -a array_identifier
declare -a array_size
declare -a array_level
runs=0;
file_res=0;

TYPE="$4";

if [[ "$TYPE" == "NM" ]]; then

    printf "Identifier\tSize1\tSize2\tSize3\tSize4\tSize7\tSize12\tSize14\tLevel1\tLevel2\tLevel3\tLevel4\tLevel7\tLevel12\tLevel14\n" >> Analysis/GeCo/$1/$3

    readarray -t array_identifier < Analysis/GeCo/$1/$2_ID.txt

    while read line; 
    do
        if [[ "$runs" == "7" ]]; then
            printf "${array_identifier[$file_res]}\t${array_size[0]}\t${array_size[1]}\t${array_size[2]}\t${array_size[3]}\t${array_size[4]}\t${array_size[5]}\t${array_size[6]}\t${array_level[0]}\t${array_level[1]}\t${array_level[2]}\t${array_level[3]}\t${array_level[4]}\t${array_level[5]}\t${array_level[6]}\n" >> Analysis/GeCo/$1/$3
            (( file_res++ ))
            runs=0;
            array_size=()
            array_level=()
        fi

        array_size[$runs]=`echo $line | cut -d ' ' -f3`
        array_level[$runs]=`echo $line | cut -d ' ' -f8`
        array_level[$runs]=`echo "scale=5; ${array_level[$runs]}/2" | bc`
        (( runs++ ))
        
    done <Analysis/GeCo/$1/geco_$2.txt

    # Input of the last set of data
    printf "${array_identifier[$file_res]}\t${array_size[0]}\t${array_size[1]}\t${array_size[2]}\t${array_size[3]}\t${array_size[4]}\t${array_size[5]}\t${array_size[6]}\t${array_level[0]}\t${array_level[1]}\t${array_level[2]}\t${array_level[3]}\t${array_level[4]}\t${array_level[5]}\t${array_level[6]}\n" >> Analysis/GeCo/$1/$3

elif [[ "$TYPE" == "PT" ]]; then

    printf "Identifier\tSize1\tSize2\tSize3\tSize4\tSize7\tSize12\tSize14\tLevel1\tLevel2\tLevel3\tLevel4\tLevel7\tLevel12\tLevel14\n" >> Analysis/AC/$1/$3

    readarray -t array_identifier < Analysis/AC/$1/$2_ID.txt

    while read line; 
    do
        if [[ "$runs" == "7" ]]; then
            printf "${array_identifier[$file_res]}\t${array_size[0]}\t${array_size[1]}\t${array_size[2]}\t${array_size[3]}\t${array_size[4]}\t${array_size[5]}\t${array_size[6]}\t${array_level[0]}\t${array_level[1]}\t${array_level[2]}\t${array_level[3]}\t${array_level[4]}\t${array_level[5]}\t${array_level[6]}\n" >> Analysis/AC/$1/$3
            (( file_res++ ))
            runs=0;
            array_size=()
            array_level=()
        fi

        array_size[$runs]=`echo $line | cut -d ' ' -f3`
        array_level[$runs]=`echo $line | cut -d ' ' -f16`
        (( runs++ ))
        
    done <Analysis/AC/$1/ac_$2.txt

    # Input of the last set of data
    printf "${array_identifier[$file_res]}\t${array_size[0]}\t${array_size[1]}\t${array_size[2]}\t${array_size[3]}\t${array_size[4]}\t${array_size[5]}\t${array_size[6]}\t${array_level[0]}\t${array_level[1]}\t${array_level[2]}\t${array_level[3]}\t${array_level[4]}\t${array_level[5]}\t${array_level[6]}\n" >> Analysis/AC/$1/$3
fi