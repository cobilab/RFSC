#!/bin/bash

rm -rf GeCo3_Output
mkdir GeCo3_Output

################## VIRUS ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Virus Genomes"

mkdir GeCo3_Output/Virus

for file in References/NCBI-Virus/*
do
    if [[ "$file" == "viral_url_donwload.txt" ]]; then
        echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
    else
        zcat $file | grep -v ">" | tr -d -c "ACGT" > TO_COMPRESS
        GeCo3 -v -l 3 TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/geco_VIRAL.txt
        rm TO_COMPRESS*
    fi
done

#for file in References/NCBI-Virus/GB_DB_VIRAL/*
#do
#    zcat $file | grep -v ">" | tr -d -c "ACGT" > TO_COMPRESS
#    GeCo3 -v -l 3 TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/geco_GB_DB_VIRAL.txt
#    rm TO_COMPRESS*
#done
#
#for file in References/NCBI-Virus/GB_DB_BACTERIA_RNA/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/geco_GB_DB_BACTERIA_RNA.txt
#    rm $file.co
#done
#
#for file in References/NCBI-Virus/GB_DB_BACTERIA_CDS/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/geco_GB_DB_BACTERIA_CDS.txt
#    rm $file.co
#done

################## BACTERIA ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Bacteria Genomes"

mkdir GeCo3_Output/Bacteria

for file in References/NCBI-Bacteria/*
do
    if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
        echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
    else
        zcat $file | grep -v ">" | tr -d -c "ACGT" > TO_COMPRESS
        GeCo3 -v -l 3 TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/geco_BACTERIA.txt
        rm TO_COMPRESS*
    fi
done

#for file in References/NCBI-Bacteria/GB_DB_BACTERIA/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/geco_GB_DB_BACTERIA.txt
#    rm $file.co
#done

#for file in References/NCBI-Bacteria/GB_DB_BACTERIA_RNA/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/geco_GB_DB_BACTERIA_RNA.txt
#    rm $file.co
#done
#
#for file in References/NCBI-Bacteria/GB_DB_BACTERIA_CDS/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/geco_GB_DB_BACTERIA_CDS.txt
#    rm $file.co
#done

################## ARCHAEA ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Archaea Genomes"

mkdir GeCo3_Output/Archaea

for file in References/NCBI-Archaea/*
do
    if [[ "$file" == "archaea_url_donwload.txt" ]]; then
        echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
    else
        zcat $file | grep -v ">" | tr -d -c "ACGT" > TO_COMPRESS
        GeCo3 -v -l 3 TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/geco_ARCHAEA.txt
        rm TO_COMPRESS*
    fi
done

#for file in References/NCBI-Archaea/GB_DB_ARCHAEA/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/geco_GB_DB_ARCHAEA.txt
#    rm $file.co
#done
#
#for file in References/NCBI-Archaea/GB_DB_ARCHAEA_CDS/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/geco_GB_DB_ARCHAEA_CDS.txt
#    rm $file.co
#done
#
#for file in References/NCBI-Archaea/GB_DB_ARCHAEA_RNA/*
#do
#    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/geco_GB_DB_ARCHAEA_RNA.txt
#    rm $file.co
#done

################## FUNGI ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Fungi Genomes"

mkdir GeCo3_Output/Fungi

for file in References/NCBI-Fungi/GB_DB_FUNGI/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/geco_GB_DB_FUNGI.txt
    rm $file.co
done

for file in References/NCBI-Fungi/GB_DB_FUNGI_CDS/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/geco_GB_DB_FUNGI_CDS.txt
    rm $file.co
done

for file in References/NCBI-Fungi/GB_DB_FUNGI_RNA/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/geco_GB_DB_FUNGI_RNA.txt
    rm $file.co
done

################## PLANT ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plant Genomes"

mkdir GeCo3_Output/Plant

for file in References/NCBI-Plant/GB_DB_PLANT/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/geco_GB_DB_PLANT.txt
    rm $file.co
done

for file in References/NCBI-Plant/GB_DB_PLANT_CDS/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/geco_GB_DB_PLANT_CDS.txt
    rm $file.co
done

for file in References/NCBI-Plant/GB_DB_PLANT_RNA/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/geco_GB_DB_PLANT_RNA.txt
    rm $file.co
done

################## PROTOZOA ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Protozoa Genomes"

mkdir GeCo3_Output/Protozoa

for file in References/NCBI-Protozoa/GB_DB_PROTOZOA/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/geco_GB_DB_PROTOZOA.txt
    rm $file.co
done

for file in References/NCBI-Protozoa/GB_DB_PROTOZOA_CDS/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/geco_GB_DB_PROTOZOA_CDS.txt
    rm $file.co
done

for file in References/NCBI-Protozoa/GB_DB_PROTOZOA_RNA/*
do
    GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/geco_GB_DB_PROTOZOA_RNA.txt
    rm $file.co
done

################## PLASTID ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plastid Genomes"

mkdir GeCo3_Output/Plastid

for file in References/NCBI-Plastid/*
do
    if [[ $file == "PLDB.fa" ]]; then
        continue
    else
        GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Plastid/geco_Plastid.txt
        rm $file.co
    fi
done

################## MITOCHONDRIAL ##################

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

mkdir GeCo3_Output/Mitochondrial

for file in References/NCBI-Mitochondrial/*
do
    if [[ $file == "MTDB.fa" ]]; then
        continue
    else
        GeCo3 -v -l 3 $file | sed '1,6d' | sed '2d' >> GeCo3_Output/Mitochondrial/geco_Mitochondrial.txt
        rm $file.co
    fi
done