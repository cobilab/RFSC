#!/bin/bash

testing_levels=(1 2 3 4 7 12 14);

################## VIRUS ##################

# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Virus Genomes"

# mkdir GeCo3_Output/Virus

# for file in References/NCBI-Virus/NM-viral/*
# do
#     if [[ "$file" == "viral_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         out_file="${out_file%".fna.gz"}"
#         printf "$out_file\n" >> GeCo3_Output/Virus/VIRAL_ID.txt

#         for l in ${testing_levels[@]}; 
#         do
#             zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Virus/TO_COMPRESS
#             GeCo3 -v -l $l GeCo3_Output/Virus/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Virus/geco_VIRAL.txt
#             rm GeCo3_Output/Virus/TO_COMPRESS*
#         done
#     fi
# done

################## BACTERIA ##################

# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Bacteria Genomes"

# mkdir GeCo3_Output/Bacteria

# for file in References/NCBI-Bacteria/*
# do
#     if [[ "$file" == "bacteria_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         zcat $file | grep -v ">" | tr -d -c "ACGT" > TO_COMPRESS
#         GeCo3 -v -l 3 TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Bacteria/geco_BACTERIA.txt
#         rm TO_COMPRESS*
#     fi
# done

# ################## ARCHAEA ##################
# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Archaea Genomes"

# mkdir GeCo3_Output/Archaea

# #for file in References/NCBI-Archaea/NM-archaea/*
# do
#     if [[ "$file" == "archaea_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         out_file="${out_file%".fna.gz"}"
#         printf "$out_file\n" >> GeCo3_Output/Archaea/ARCHAEA_ID.txt

#         for l in ${testing_levels[@]}; 
#         do
#             zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Archaea/TO_COMPRESS
#             GeCo3 -v -l $l GeCo3_Output/Archaea/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Archaea/geco_ARCHAEA.txt
#             rm GeCo3_Output/Archaea/TO_COMPRESS*
#         done
#     fi
# done

# ################## FUNGI ##################
# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Fungi Genomes"

# mkdir GeCo3_Output/Fungi

# for file in References/NCBI-Fungi/NM-fungi/*
# do
#     if [[ "$file" == "fungi_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         out_file="${out_file%".fna.gz"}"
#         printf "$out_file\n" >> GeCo3_Output/Fungi/FUNGI_ID.txt

#         for l in ${testing_levels[@]}; 
#         do
#             zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Fungi/TO_COMPRESS
#             GeCo3 -v -l $l GeCo3_Output/Fungi/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Fungi/geco_FUNGI.txt
#             rm GeCo3_Output/Fungi/TO_COMPRESS*
#         done
#     fi
# done

# ################## PLANT ##################
# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plant Genomes"

# mkdir GeCo3_Output/Plant

# for file in References/NCBI-Plant/NM-plant/*
# do
#     if [[ "$file" == "plant_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         out_file="${out_file%".fna.gz"}"
#         printf "$out_file\n" >> GeCo3_Output/Plant/PLANT_ID.txt

#         for l in ${testing_levels[@]}; 
#         do
#             zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plant/TO_COMPRESS
#             GeCo3 -v -l $l GeCo3_Output/Plant/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plant/geco_PLANT.txt
#             rm GeCo3_Output/Plant/TO_COMPRESS*
#         done
#     fi
# done

# ################## PROTOZOA ##################
# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Protozoa Genomes"

# mkdir GeCo3_Output/Protozoa

# for file in References/NCBI-Protozoa/NM-protozoa/*
# do
#     if [[ "$file" == "protozoa_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         out_file="${out_file%".fna.gz"}"
#         printf "$out_file\n" >> GeCo3_Output/Protozoa/PROTOZOA_ID.txt

#         for l in ${testing_levels[@]}; 
#         do
#             zcat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Protozoa/TO_COMPRESS
#             GeCo3 -v -l $l GeCo3_Output/Protozoa/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Protozoa/geco_PROTOZOA.txt
#             rm GeCo3_Output/Protozoa/TO_COMPRESS*
#         done
#     fi
# done

################## PLASTID ##################
# echo -e "\033[1;34m[RFSC]\033[0m Start parsing Plastid Genomes"
# path="/media/alexloure/T7Touch/NCBI-Plastid"
# gunzip -k $path/DB-plastid.fa.gz
# mkdir $path/NM-plastid
# awk '/>/{filename="/media/alexloure/T7Touch/NCBI-Plastid/NM-plastid/"NR".fna"}; {print >filename; close(filename)}' $path/DB-plastid.fa

# echo -e "\033[1;34m[RFSC]\033[0m Start compressing Plastid Genomes"

# mkdir GeCo3_Output/Plastid

# for file in References/NCBI-Plastid/NM-plastid/*
# do
#     if [[ "$file" == "plastid_url_donwload.txt" ]]; then
#         echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
#     else
#         out_file=$(basename $file);
#         out_file="${out_file%".fna"}"
#         printf "$out_file\n" >> GeCo3_Output/Plastid/PLASTID_ID.txt

#         for l in ${testing_levels[@]}; 
#         do
#             cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Plastid/TO_COMPRESS
#             GeCo3 -v -l $l GeCo3_Output/Plastid/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Plastid/geco_PLASTID.txt
#             rm GeCo3_Output/Plastid/TO_COMPRESS*
#         done
#     fi
# done

# ################## MITOCHONDRIAL ##################
echo -e "\033[1;34m[RFSC]\033[0m Start parsing Mitochondrial Genomes"
path="/media/alexloure/T7Touch/NCBI-Mitochondrial"
gunzip -k $path/DB-mitochondrion.fa.gz
mkdir $path/NM-mitochondrion
awk '/>/{filename="/media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/"NR".fna"}; {print >filename; close(filename)}' $path/DB-mitochondrion.fa

echo -e "\033[1;34m[RFSC]\033[0m Start compressing Mitochondrial Genomes"

mkdir GeCo3_Output/Mitochondrial

#for file in References/NCBI-Mitochondrial/NM-mitochondrion/*
for file in /media/alexloure/T7Touch/NCBI-Mitochondrial/NM-mitochondrion/*
do
    if [[ "$file" == "mitochondrion_url_donwload.txt" ]]; then
        echo -e "\033[1;34m[RFSC]\033[0m Jumping index file ($file)";
    else
        out_file=$(basename $file);
        out_file="${out_file%".fna"}"
        printf "$out_file\n" >> GeCo3_Output/Mitochondrial/MITOCHONDRIAL_ID.txt

        for l in ${testing_levels[@]}; 
        do
            cat $file | grep -v ">" | tr -d -c "ACGT" > GeCo3_Output/Mitochondrial/TO_COMPRESS
            GeCo3 -v -l $l GeCo3_Output/Mitochondrial/TO_COMPRESS | sed '1,6d' | sed '2d' >> GeCo3_Output/Mitochondrial/geco_MITOCHONDRIAL.txt
            rm GeCo3_Output/Mitochondrial/TO_COMPRESS*
        done
    fi
done