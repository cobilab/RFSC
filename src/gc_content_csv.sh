#!/bin/bash

# run: ./src/gc_content_csv.sh Virus VIRAL gc_content_Viral.csv
#      ./src/gc_content_csv.sh Bacteria BACTERIA gc_content_Bacteria.csv
#      ./src/gc_content_csv.sh Fungi FUNGI gc_content_Fungi.csv
#      ./src/gc_content_csv.sh Archaea ARCHAEA gc_content_Archaea.csv
#      ./src/gc_content_csv.sh Plant PLANT gc_content_Plant.csv
#      ./src/gc_content_csv.sh Protozoa PROTOZOA gc_content_Protozoa.csv
#      ./src/gc_content_csv.sh Plastid PLASTID gc_content_Plastid.csv
#      ./src/gc_content_csv.sh Mitochondrial MITOCHONDRIAL gc_content_Mitochondrial.csv

declare -a array_identifier
GC_VALUE="";
ATGC_VALUE="";
PERCENTAGE=""; # Between [0..1]

file_res=0;

printf "Identifier\tGC\tATGC\tPercentage\n" >> Analysis/GCcontent/$1/$3

readarray -t array_identifier < Analysis/GCcontent/$1/$2_ID.txt

while read line; 
do
    GC_VALUE=`echo $line | cut -d ' ' -f2`
    ATGC_VALUE=`echo $line | cut -d ' ' -f4`
    PERCENTAGE=`echo $line | cut -d ' ' -f6`

    printf "${array_identifier[$file_res]}\t$GC_VALUE\t$ATGC_VALUE\t$PERCENTAGE\n" >> Analysis/GCcontent/$1/$3
    (( file_res++ ))
    
done <Analysis/GCcontent/$1/GC_content_$2.txt