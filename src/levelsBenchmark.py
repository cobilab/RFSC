from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv

## Virus Values
virus_nm_size1, virus_nm_size2, virus_nm_size3, virus_nm_size4, virus_nm_size7, virus_nm_size12, virus_nm_size14 = [], [], [], [], [], [], []
virus_pt_size1, virus_pt_size2, virus_pt_size3, virus_pt_size4, virus_pt_size7, virus_pt_size12, virus_pt_size14 = [], [], [], [], [], [], []

## Bacterias Values
bacterias_nm_size1, bacterias_nm_size2, bacterias_nm_size3, bacterias_nm_size4, bacterias_nm_size7, bacterias_nm_size12, bacterias_nm_size14 = [], [], [], [], [], [], []
bacterias_pt_size1, bacterias_pt_size2, bacterias_pt_size3, bacterias_pt_size4, bacterias_pt_size7, bacterias_pt_size12, bacterias_pt_size14 = [], [], [], [], [], [], []

## Archaeas Values
archaeas_nm_size1, archaeas_nm_size2, archaeas_nm_size3, archaeas_nm_size4, archaeas_nm_size7, archaeas_nm_size12, archaeas_nm_size14 = [], [], [], [], [], [], []
archaeas_pt_size1, archaeas_pt_size2, archaeas_pt_size3, archaeas_pt_size4, archaeas_pt_size7, archaeas_pt_size12, archaeas_pt_size14 = [], [], [], [], [], [], []

## Fungi Values
fungi_nm_size1, fungi_nm_size2, fungi_nm_size3, fungi_nm_size4, fungi_nm_size7, fungi_nm_size12, fungi_nm_size14 = [], [], [], [], [], [], []
fungi_pt_size1, fungi_pt_size2, fungi_pt_size3, fungi_pt_size4, fungi_pt_size7, fungi_pt_size12, fungi_pt_size14 = [], [], [], [], [], [], []

## Plant Values
plant_nm_size1, plant_nm_size2, plant_nm_size3, plant_nm_size4, plant_nm_size7, plant_nm_size12, plant_nm_size14 = [], [], [], [], [], [], []
plant_pt_size1, plant_pt_size2, plant_pt_size3, plant_pt_size4, plant_pt_size7, plant_pt_size12, plant_pt_size14 = [], [], [], [], [], [], []

## Protozoa Values
protozoa_nm_size1, protozoa_nm_size2, protozoa_nm_size3, protozoa_nm_size4, protozoa_nm_size7, protozoa_nm_size12, protozoa_nm_size14 = [], [], [], [], [], [], []
protozoa_pt_size1, protozoa_pt_size2, protozoa_pt_size3, protozoa_pt_size4, protozoa_pt_size7, protozoa_pt_size12, protozoa_pt_size14 = [], [], [], [], [], [], []

## Mitochondrial Values
mito_nm_size1, mito_nm_size2, mito_nm_size3, mito_nm_size4, mito_nm_size7, mito_nm_size12, mito_nm_size14 = [], [], [], [], [], [], []
mito_pt_size1, mito_pt_size2, mito_pt_size3, mito_pt_size4, mito_pt_size7, mito_pt_size12, mito_pt_size14 = [], [], [], [], [], [], []

## Plastid Values
plastid_nm_size1, plastid_nm_size2, plastid_nm_size3, plastid_nm_size4, plastid_nm_size7, plastid_nm_size12, plastid_nm_size14 = [], [], [], [], [], [], []
plastid_pt_size1, plastid_pt_size2, plastid_pt_size3, plastid_pt_size4, plastid_pt_size7, plastid_pt_size12, plastid_pt_size14 = [], [], [], [], [], [], []

# Virus CSV files
with open('GeCo3_Output/Virus/NM/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_nm_size1.append(row[0].split("\t")[1])
        virus_nm_size2.append(row[0].split("\t")[2])
        virus_nm_size3.append(row[0].split("\t")[3])
        virus_nm_size4.append(row[0].split("\t")[4])
        virus_nm_size7.append(row[0].split("\t")[5])
        virus_nm_size12.append(row[0].split("\t")[6])
        virus_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Virus/PT/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_pt_size1.append(row[0].split("\t")[1])
        virus_pt_size2.append(row[0].split("\t")[2])
        virus_pt_size3.append(row[0].split("\t")[3])
        virus_pt_size4.append(row[0].split("\t")[4])
        virus_pt_size7.append(row[0].split("\t")[5])
        virus_pt_size12.append(row[0].split("\t")[6])
        virus_pt_size14.append(row[0].split("\t")[7])

'''
# Bacterias CSV files
with open('GeCo3_Output/Bacteria/NM/geco3_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacterias_nm_size1.append(row[0].split("\t")[1])
        bacterias_nm_size2.append(row[0].split("\t")[2])
        bacterias_nm_size3.append(row[0].split("\t")[3])
        bacterias_nm_size4.append(row[0].split("\t")[4])
        bacterias_nm_size7.append(row[0].split("\t")[5])
        bacterias_nm_size12.append(row[0].split("\t")[6])
        bacterias_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Bacteria/PT/geco3_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacterias_pt_size1.append(row[0].split("\t")[1])
        bacterias_pt_size2.append(row[0].split("\t")[2])
        bacterias_pt_size3.append(row[0].split("\t")[3])
        bacterias_pt_size4.append(row[0].split("\t")[4])
        bacterias_pt_size7.append(row[0].split("\t")[5])
        bacterias_pt_size12.append(row[0].split("\t")[6])
        bacterias_pt_size14.append(row[0].split("\t")[7])
'''

# Archaeas CSV files
with open('GeCo3_Output/Archaea/NM/geco3_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaeas_nm_size1.append(row[0].split("\t")[1])
        archaeas_nm_size2.append(row[0].split("\t")[2])
        archaeas_nm_size3.append(row[0].split("\t")[3])
        archaeas_nm_size4.append(row[0].split("\t")[4])
        archaeas_nm_size7.append(row[0].split("\t")[5])
        archaeas_nm_size12.append(row[0].split("\t")[6])
        archaeas_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Archaea/PT/geco3_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaeas_pt_size1.append(row[0].split("\t")[1])
        archaeas_pt_size2.append(row[0].split("\t")[2])
        archaeas_pt_size3.append(row[0].split("\t")[3])
        archaeas_pt_size4.append(row[0].split("\t")[4])
        archaeas_pt_size7.append(row[0].split("\t")[5])
        archaeas_pt_size12.append(row[0].split("\t")[6])
        archaeas_pt_size14.append(row[0].split("\t")[7])

# Fungis CSV files
with open('GeCo3_Output/Fungi/NM/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_nm_size1.append(row[0].split("\t")[1])
        fungi_nm_size2.append(row[0].split("\t")[2])
        fungi_nm_size3.append(row[0].split("\t")[3])
        fungi_nm_size4.append(row[0].split("\t")[4])
        fungi_nm_size7.append(row[0].split("\t")[5])
        fungi_nm_size12.append(row[0].split("\t")[6])
        fungi_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Fungi/PT/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_pt_size1.append(row[0].split("\t")[1])
        fungi_pt_size2.append(row[0].split("\t")[2])
        fungi_pt_size3.append(row[0].split("\t")[3])
        fungi_pt_size4.append(row[0].split("\t")[4])
        fungi_pt_size7.append(row[0].split("\t")[5])
        fungi_pt_size12.append(row[0].split("\t")[6])
        fungi_pt_size14.append(row[0].split("\t")[7])

# Plants CSV files
with open('GeCo3_Output/Plant/NM/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_nm_size1.append(row[0].split("\t")[1])
        plant_nm_size2.append(row[0].split("\t")[2])
        plant_nm_size3.append(row[0].split("\t")[3])
        plant_nm_size4.append(row[0].split("\t")[4])
        plant_nm_size7.append(row[0].split("\t")[5])
        plant_nm_size12.append(row[0].split("\t")[6])
        plant_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Plant/PT/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_pt_size1.append(row[0].split("\t")[1])
        plant_pt_size2.append(row[0].split("\t")[2])
        plant_pt_size3.append(row[0].split("\t")[3])
        plant_pt_size4.append(row[0].split("\t")[4])
        plant_pt_size7.append(row[0].split("\t")[5])
        plant_pt_size12.append(row[0].split("\t")[6])
        plant_pt_size14.append(row[0].split("\t")[7])

# Protozoa CSV files
with open('GeCo3_Output/Protozoa/NM/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_nm_size1.append(row[0].split("\t")[1])
        protozoa_nm_size2.append(row[0].split("\t")[2])
        protozoa_nm_size3.append(row[0].split("\t")[3])
        protozoa_nm_size4.append(row[0].split("\t")[4])
        protozoa_nm_size7.append(row[0].split("\t")[5])
        protozoa_nm_size12.append(row[0].split("\t")[6])
        protozoa_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Protozoa/PT/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_pt_size1.append(row[0].split("\t")[1])
        protozoa_pt_size2.append(row[0].split("\t")[2])
        protozoa_pt_size3.append(row[0].split("\t")[3])
        protozoa_pt_size4.append(row[0].split("\t")[4])
        protozoa_pt_size7.append(row[0].split("\t")[5])
        protozoa_pt_size12.append(row[0].split("\t")[6])
        protozoa_pt_size14.append(row[0].split("\t")[7])

# Mitochondrial CSV files
with open('GeCo3_Output/Mitochondrial/NM/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_nm_size1.append(row[0].split("\t")[1])
        mito_nm_size2.append(row[0].split("\t")[2])
        mito_nm_size3.append(row[0].split("\t")[3])
        mito_nm_size4.append(row[0].split("\t")[4])
        mito_nm_size7.append(row[0].split("\t")[5])
        mito_nm_size12.append(row[0].split("\t")[6])
        mito_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Mitochondrial/PT/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_pt_size1.append(row[0].split("\t")[1])
        mito_pt_size2.append(row[0].split("\t")[2])
        mito_pt_size3.append(row[0].split("\t")[3])
        mito_pt_size4.append(row[0].split("\t")[4])
        mito_pt_size7.append(row[0].split("\t")[5])
        mito_pt_size12.append(row[0].split("\t")[6])
        mito_pt_size14.append(row[0].split("\t")[7])

# Plastid CSV files
with open('GeCo3_Output/Plastid/NM/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_nm_size1.append(row[0].split("\t")[1])
        plastid_nm_size2.append(row[0].split("\t")[2])
        plastid_nm_size3.append(row[0].split("\t")[3])
        plastid_nm_size4.append(row[0].split("\t")[4])
        plastid_nm_size7.append(row[0].split("\t")[5])
        plastid_nm_size12.append(row[0].split("\t")[6])
        plastid_nm_size14.append(row[0].split("\t")[7])

with open('GeCo3_Output/Plastid/PT/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_pt_size1.append(row[0].split("\t")[1])
        plastid_pt_size2.append(row[0].split("\t")[2])
        plastid_pt_size3.append(row[0].split("\t")[3])
        plastid_pt_size4.append(row[0].split("\t")[4])
        plastid_pt_size7.append(row[0].split("\t")[5])
        plastid_pt_size12.append(row[0].split("\t")[6])
        plastid_pt_size14.append(row[0].split("\t")[7])

'''
    Convert sizes from string to int
'''

## Virus conversion
virus_nm_size1 = np.array(virus_nm_size1[1:], dtype=int)
virus_nm_size2 = np.array(virus_nm_size2[1:], dtype=int)
virus_nm_size3 = np.array(virus_nm_size3[1:], dtype=int)
virus_nm_size4 = np.array(virus_nm_size4[1:], dtype=int)
virus_nm_size7 = np.array(virus_nm_size7[1:], dtype=int)
virus_nm_size12 = np.array(virus_nm_size12[1:], dtype=int)
virus_nm_size14 = np.array(virus_nm_size14[1:], dtype=int)

virus_pt_size1 = np.array(virus_pt_size1[1:], dtype=int)
virus_pt_size2 = np.array(virus_pt_size2[1:], dtype=int)
virus_pt_size3 = np.array(virus_pt_size3[1:], dtype=int)
virus_pt_size4 = np.array(virus_pt_size4[1:], dtype=int)
virus_pt_size7 = np.array(virus_pt_size7[1:], dtype=int)
virus_pt_size12 = np.array(virus_pt_size12[1:], dtype=int)
virus_pt_size14 = np.array(virus_pt_size14[1:], dtype=int)

## Bacterias conversion
bacterias_nm_size1 = np.array(bacterias_nm_size1[1:], dtype=int)
bacterias_nm_size2 = np.array(bacterias_nm_size2[1:], dtype=int)
bacterias_nm_size3 = np.array(bacterias_nm_size3[1:], dtype=int)
bacterias_nm_size4 = np.array(bacterias_nm_size4[1:], dtype=int)
bacterias_nm_size7 = np.array(bacterias_nm_size7[1:], dtype=int)
bacterias_nm_size12 = np.array(bacterias_nm_size12[1:], dtype=int)
bacterias_nm_size14 = np.array(bacterias_nm_size14[1:], dtype=int)

bacterias_pt_size1 = np.array(bacterias_pt_size1[1:], dtype=int)
bacterias_pt_size2 = np.array(bacterias_pt_size2[1:], dtype=int)
bacterias_pt_size3 = np.array(bacterias_pt_size3[1:], dtype=int)
bacterias_pt_size4 = np.array(bacterias_pt_size4[1:], dtype=int)
bacterias_pt_size7 = np.array(bacterias_pt_size7[1:], dtype=int)
bacterias_pt_size12 = np.array(bacterias_pt_size12[1:], dtype=int)
bacterias_pt_size14 = np.array(bacterias_pt_size14[1:], dtype=int)

## Archaeas conversion
archaeas_nm_size1 = np.array(archaeas_nm_size1[1:], dtype=int)
archaeas_nm_size2 = np.array(archaeas_nm_size2[1:], dtype=int)
archaeas_nm_size3 = np.array(archaeas_nm_size3[1:], dtype=int)
archaeas_nm_size4 = np.array(archaeas_nm_size4[1:], dtype=int)
archaeas_nm_size7 = np.array(archaeas_nm_size7[1:], dtype=int)
archaeas_nm_size12 = np.array(archaeas_nm_size12[1:], dtype=int)
archaeas_nm_size14 = np.array(archaeas_nm_size14[1:], dtype=int)

archaeas_pt_size1 = np.array(archaeas_pt_size1[1:], dtype=int)
archaeas_pt_size2 = np.array(archaeas_pt_size2[1:], dtype=int)
archaeas_pt_size3 = np.array(archaeas_pt_size3[1:], dtype=int)
archaeas_pt_size4 = np.array(archaeas_pt_size4[1:], dtype=int)
archaeas_pt_size7 = np.array(archaeas_pt_size7[1:], dtype=int)
archaeas_pt_size12 = np.array(archaeas_pt_size12[1:], dtype=int)
archaeas_pt_size14 = np.array(archaeas_pt_size14[1:], dtype=int)

## Fungis conversion
fungi_nm_size1 = np.array(fungi_nm_size1[1:], dtype=int)
fungi_nm_size2 = np.array(fungi_nm_size2[1:], dtype=int)
fungi_nm_size3 = np.array(fungi_nm_size3[1:], dtype=int)
fungi_nm_size4 = np.array(fungi_nm_size4[1:], dtype=int)
fungi_nm_size7 = np.array(fungi_nm_size7[1:], dtype=int)
fungi_nm_size12 = np.array(fungi_nm_size12[1:], dtype=int)
fungi_nm_size14 = np.array(fungi_nm_size14[1:], dtype=int)

fungi_pt_size1 = np.array(fungi_pt_size1[1:], dtype=int)
fungi_pt_size2 = np.array(fungi_pt_size2[1:], dtype=int)
fungi_pt_size3 = np.array(fungi_pt_size3[1:], dtype=int)
fungi_pt_size4 = np.array(fungi_pt_size4[1:], dtype=int)
fungi_pt_size7 = np.array(fungi_pt_size7[1:], dtype=int)
fungi_pt_size12 = np.array(fungi_pt_size12[1:], dtype=int)
fungi_pt_size14 = np.array(fungi_pt_size14[1:], dtype=int)

## Plant conversion
plant_nm_size1 = np.array(plant_nm_size1[1:], dtype=int)
plant_nm_size2 = np.array(plant_nm_size2[1:], dtype=int)
plant_nm_size3 = np.array(plant_nm_size3[1:], dtype=int)
plant_nm_size4 = np.array(plant_nm_size4[1:], dtype=int)
plant_nm_size7 = np.array(plant_nm_size7[1:], dtype=int)
plant_nm_size12 = np.array(plant_nm_size12[1:], dtype=int)
plant_nm_size14 = np.array(plant_nm_size14[1:], dtype=int)

plant_pt_size1 = np.array(plant_pt_size1[1:], dtype=int)
plant_pt_size2 = np.array(plant_pt_size2[1:], dtype=int)
plant_pt_size3 = np.array(plant_pt_size3[1:], dtype=int)
plant_pt_size4 = np.array(plant_pt_size4[1:], dtype=int)
plant_pt_size7 = np.array(plant_pt_size7[1:], dtype=int)
plant_pt_size12 = np.array(plant_pt_size12[1:], dtype=int)
plant_pt_size14 = np.array(plant_pt_size14[1:], dtype=int)

## Protozoa conversion
protozoa_nm_size1 = np.array(protozoa_nm_size1[1:], dtype=int)
protozoa_nm_size2 = np.array(protozoa_nm_size2[1:], dtype=int)
protozoa_nm_size3 = np.array(protozoa_nm_size3[1:], dtype=int)
protozoa_nm_size4 = np.array(protozoa_nm_size4[1:], dtype=int)
protozoa_nm_size7 = np.array(protozoa_nm_size7[1:], dtype=int)
protozoa_nm_size12 = np.array(protozoa_nm_size12[1:], dtype=int)
protozoa_nm_size14 = np.array(protozoa_nm_size14[1:], dtype=int)

protozoa_pt_size1 = np.array(protozoa_pt_size1[1:], dtype=int)
protozoa_pt_size2 = np.array(protozoa_pt_size2[1:], dtype=int)
protozoa_pt_size3 = np.array(protozoa_pt_size3[1:], dtype=int)
protozoa_pt_size4 = np.array(protozoa_pt_size4[1:], dtype=int)
protozoa_pt_size7 = np.array(protozoa_pt_size7[1:], dtype=int)
protozoa_pt_size12 = np.array(protozoa_pt_size12[1:], dtype=int)
protozoa_pt_size14 = np.array(protozoa_pt_size14[1:], dtype=int)

## Mitochondrial conversion
mito_nm_size1 = np.array(mito_nm_size1[1:], dtype=int)
mito_nm_size2 = np.array(mito_nm_size2[1:], dtype=int)
mito_nm_size3 = np.array(mito_nm_size3[1:], dtype=int)
mito_nm_size4 = np.array(mito_nm_size4[1:], dtype=int)
mito_nm_size7 = np.array(mito_nm_size7[1:], dtype=int)
mito_nm_size12 = np.array(mito_nm_size12[1:], dtype=int)
mito_nm_size14 = np.array(mito_nm_size14[1:], dtype=int)

mito_pt_size1 = np.array(mito_pt_size1[1:], dtype=int)
mito_pt_size2 = np.array(mito_pt_size2[1:], dtype=int)
mito_pt_size3 = np.array(mito_pt_size3[1:], dtype=int)
mito_pt_size4 = np.array(mito_pt_size4[1:], dtype=int)
mito_pt_size7 = np.array(mito_pt_size7[1:], dtype=int)
mito_pt_size12 = np.array(mito_pt_size12[1:], dtype=int)
mito_pt_size14 = np.array(mito_pt_size14[1:], dtype=int)

## Plastid conversion
plastid_nm_size1 = np.array(plastid_nm_size1[1:], dtype=int)
plastid_nm_size2 = np.array(plastid_nm_size2[1:], dtype=int)
plastid_nm_size3 = np.array(plastid_nm_size3[1:], dtype=int)
plastid_nm_size4 = np.array(plastid_nm_size4[1:], dtype=int)
plastid_nm_size7 = np.array(plastid_nm_size7[1:], dtype=int)
plastid_nm_size12 = np.array(plastid_nm_size12[1:], dtype=int)
plastid_nm_size14 = np.array(plastid_nm_size14[1:], dtype=int)

plastid_pt_size1 = np.array(plastid_pt_size1[1:], dtype=int)
plastid_pt_size2 = np.array(plastid_pt_size2[1:], dtype=int)
plastid_pt_size3 = np.array(plastid_pt_size3[1:], dtype=int)
plastid_pt_size4 = np.array(plastid_pt_size4[1:], dtype=int)
plastid_pt_size7 = np.array(plastid_pt_size7[1:], dtype=int)
plastid_pt_size12 = np.array(plastid_pt_size12[1:], dtype=int)
plastid_pt_size14 = np.array(plastid_pt_size14[1:], dtype=int)

'''
    Calculate cumulative sizes for every level of each domain
'''

# sum(Virus)
virus_nm_level1_cumulative = sum(virus_nm_size1)
virus_nm_level2_cumulative = sum(virus_nm_size2)
virus_nm_level3_cumulative = sum(virus_nm_size3)
virus_nm_level4_cumulative = sum(virus_nm_size4)
virus_nm_level7_cumulative = sum(virus_nm_size7)
virus_nm_level12_cumulative = sum(virus_nm_size12)
virus_nm_level14_cumulative = sum(virus_nm_size14)

virus_pt_level1_cumulative = sum(virus_pt_size1)
virus_pt_level2_cumulative = sum(virus_pt_size2)
virus_pt_level3_cumulative = sum(virus_pt_size3)
virus_pt_level4_cumulative = sum(virus_pt_size4)
virus_pt_level7_cumulative = sum(virus_pt_size7)
virus_pt_level12_cumulative = sum(virus_pt_size12)
virus_pt_level14_cumulative = sum(virus_pt_size14)

# sum(Bacteria)
bacterias_nm_level1_cumulative = sum(bacterias_nm_size1)
bacterias_nm_level2_cumulative = sum(bacterias_nm_size2)
bacterias_nm_level3_cumulative = sum(bacterias_nm_size3)
bacterias_nm_level4_cumulative = sum(bacterias_nm_size4)
bacterias_nm_level7_cumulative = sum(bacterias_nm_size7)
bacterias_nm_level12_cumulative = sum(bacterias_nm_size12)
bacterias_nm_level14_cumulative = sum(bacterias_nm_size14)

bacterias_pt_level1_cumulative = sum(bacterias_pt_size1)
bacterias_pt_level2_cumulative = sum(bacterias_pt_size2)
bacterias_pt_level3_cumulative = sum(bacterias_pt_size3)
bacterias_pt_level4_cumulative = sum(bacterias_pt_size4)
bacterias_pt_level7_cumulative = sum(bacterias_pt_size7)
bacterias_pt_level12_cumulative = sum(bacterias_pt_size12)
bacterias_pt_level14_cumulative = sum(bacterias_pt_size14)

# sum(Archaea)
archaeas_nm_level1_cumulative = sum(archaeas_nm_size1)
archaeas_nm_level2_cumulative = sum(archaeas_nm_size2)
archaeas_nm_level3_cumulative = sum(archaeas_nm_size3)
archaeas_nm_level4_cumulative = sum(archaeas_nm_size4)
archaeas_nm_level7_cumulative = sum(archaeas_nm_size7)
archaeas_nm_level12_cumulative = sum(archaeas_nm_size12)
archaeas_nm_level14_cumulative = sum(archaeas_nm_size14)

archaeas_pt_level1_cumulative = sum(archaeas_pt_size1)
archaeas_pt_level2_cumulative = sum(archaeas_pt_size2)
archaeas_pt_level3_cumulative = sum(archaeas_pt_size3)
archaeas_pt_level4_cumulative = sum(archaeas_pt_size4)
archaeas_pt_level7_cumulative = sum(archaeas_pt_size7)
archaeas_pt_level12_cumulative = sum(archaeas_pt_size12)
archaeas_pt_level14_cumulative = sum(archaeas_pt_size14)

# sum(Fungi)
fungi_nm_level1_cumulative = sum(fungi_nm_size1)
fungi_nm_level2_cumulative = sum(fungi_nm_size2)
fungi_nm_level3_cumulative = sum(fungi_nm_size3)
fungi_nm_level4_cumulative = sum(fungi_nm_size4)
fungi_nm_level7_cumulative = sum(fungi_nm_size7)
fungi_nm_level12_cumulative = sum(fungi_nm_size12)
fungi_nm_level14_cumulative = sum(fungi_nm_size14)

fungi_pt_level1_cumulative = sum(fungi_pt_size1)
fungi_pt_level2_cumulative = sum(fungi_pt_size2)
fungi_pt_level3_cumulative = sum(fungi_pt_size3)
fungi_pt_level4_cumulative = sum(fungi_pt_size4)
fungi_pt_level7_cumulative = sum(fungi_pt_size7)
fungi_pt_level12_cumulative = sum(fungi_pt_size12)
fungi_pt_level14_cumulative = sum(fungi_pt_size14)

# sum(Plant)
plant_nm_level1_cumulative = sum(plant_nm_size1)
plant_nm_level2_cumulative = sum(plant_nm_size2)
plant_nm_level3_cumulative = sum(plant_nm_size3)
plant_nm_level4_cumulative = sum(plant_nm_size4)
plant_nm_level7_cumulative = sum(plant_nm_size7)
plant_nm_level12_cumulative = sum(plant_nm_size12)
plant_nm_level14_cumulative = sum(plant_nm_size14)

plant_pt_level1_cumulative = sum(plant_pt_size1)
plant_pt_level2_cumulative = sum(plant_pt_size2)
plant_pt_level3_cumulative = sum(plant_pt_size3)
plant_pt_level4_cumulative = sum(plant_pt_size4)
plant_pt_level7_cumulative = sum(plant_pt_size7)
plant_pt_level12_cumulative = sum(plant_pt_size12)
plant_pt_level14_cumulative = sum(plant_pt_size14)

# sum(Protozoa)
protozoa_nm_level1_cumulative = sum(protozoa_nm_size1)
protozoa_nm_level2_cumulative = sum(protozoa_nm_size2)
protozoa_nm_level3_cumulative = sum(protozoa_nm_size3)
protozoa_nm_level4_cumulative = sum(protozoa_nm_size4)
protozoa_nm_level7_cumulative = sum(protozoa_nm_size7)
protozoa_nm_level12_cumulative = sum(protozoa_nm_size12)
protozoa_nm_level14_cumulative = sum(protozoa_nm_size14)

protozoa_pt_level1_cumulative = sum(protozoa_pt_size1)
protozoa_pt_level2_cumulative = sum(protozoa_pt_size2)
protozoa_pt_level3_cumulative = sum(protozoa_pt_size3)
protozoa_pt_level4_cumulative = sum(protozoa_pt_size4)
protozoa_pt_level7_cumulative = sum(protozoa_pt_size7)
protozoa_pt_level12_cumulative = sum(protozoa_pt_size12)
protozoa_pt_level14_cumulative = sum(protozoa_pt_size14)

# sum(Mitochondrial)
mito_nm_level1_cumulative = sum(mito_nm_size1)
mito_nm_level2_cumulative = sum(mito_nm_size2)
mito_nm_level3_cumulative = sum(mito_nm_size3)
mito_nm_level4_cumulative = sum(mito_nm_size4)
mito_nm_level7_cumulative = sum(mito_nm_size7)
mito_nm_level12_cumulative = sum(mito_nm_size12)
mito_nm_level14_cumulative = sum(mito_nm_size14)

mito_pt_level1_cumulative = sum(mito_pt_size1)
mito_pt_level2_cumulative = sum(mito_pt_size2)
mito_pt_level3_cumulative = sum(mito_pt_size3)
mito_pt_level4_cumulative = sum(mito_pt_size4)
mito_pt_level7_cumulative = sum(mito_pt_size7)
mito_pt_level12_cumulative = sum(mito_pt_size12)
mito_pt_level14_cumulative = sum(mito_pt_size14)

# sum(Plastid)
plastid_nm_level1_cumulative = sum(plastid_nm_size1)
plastid_nm_level2_cumulative = sum(plastid_nm_size2)
plastid_nm_level3_cumulative = sum(plastid_nm_size3)
plastid_nm_level4_cumulative = sum(plastid_nm_size4)
plastid_nm_level7_cumulative = sum(plastid_nm_size7)
plastid_nm_level12_cumulative = sum(plastid_nm_size12)
plastid_nm_level14_cumulative = sum(plastid_nm_size14)

plastid_pt_level1_cumulative = sum(plastid_pt_size1)
plastid_pt_level2_cumulative = sum(plastid_pt_size2)
plastid_pt_level3_cumulative = sum(plastid_pt_size3)
plastid_pt_level4_cumulative = sum(plastid_pt_size4)
plastid_pt_level7_cumulative = sum(plastid_pt_size7)
plastid_pt_level12_cumulative = sum(plastid_pt_size12)
plastid_pt_level14_cumulative = sum(plastid_pt_size14)


X = ['1', '2', '3', '4', '7', '12', '14']
y_virus_nm = (virus_nm_level1_cumulative, virus_nm_level2_cumulative, virus_nm_level3_cumulative, virus_nm_level4_cumulative, virus_nm_level7_cumulative, virus_nm_level12_cumulative, virus_nm_level14_cumulative)
y_bacterias_nm = (bacterias_nm_level1_cumulative, bacterias_nm_level2_cumulative, bacterias_nm_level3_cumulative, bacterias_nm_level4_cumulative, bacterias_nm_level7_cumulative, bacterias_nm_level12_cumulative, bacterias_nm_level14_cumulative)
y_archaeas_nm = (archaeas_nm_level1_cumulative, archaeas_nm_level2_cumulative, archaeas_nm_level3_cumulative, archaeas_nm_level4_cumulative, archaeas_nm_level7_cumulative, archaeas_nm_level12_cumulative, archaeas_nm_level14_cumulative)
y_fungi_nm = (fungi_nm_level1_cumulative, fungi_nm_level2_cumulative, fungi_nm_level3_cumulative, fungi_nm_level4_cumulative, fungi_nm_level7_cumulative, fungi_nm_level12_cumulative, fungi_nm_level14_cumulative)
y_plant_nm = (plant_nm_level1_cumulative, plant_nm_level2_cumulative, plant_nm_level3_cumulative, plant_nm_level4_cumulative, plant_nm_level7_cumulative, plant_nm_level12_cumulative, plant_nm_level14_cumulative)
y_protozoa_nm = (protozoa_nm_level1_cumulative, protozoa_nm_level2_cumulative, protozoa_nm_level3_cumulative, protozoa_nm_level4_cumulative, protozoa_nm_level7_cumulative, protozoa_nm_level12_cumulative, protozoa_nm_level14_cumulative)
y_mito_nm = (mito_nm_level1_cumulative, mito_nm_level2_cumulative, mito_nm_level3_cumulative, mito_nm_level4_cumulative, mito_nm_level7_cumulative, mito_nm_level12_cumulative, mito_nm_level14_cumulative)
y_plastid_nm = (plastid_nm_level1_cumulative, plastid_nm_level2_cumulative, plastid_nm_level3_cumulative, plastid_nm_level4_cumulative, plastid_nm_level7_cumulative, plastid_nm_level12_cumulative, plastid_nm_level14_cumulative)

y_virus_pt = (virus_pt_level1_cumulative, virus_pt_level2_cumulative, virus_pt_level3_cumulative, virus_pt_level4_cumulative, virus_pt_level7_cumulative, virus_pt_level12_cumulative, virus_pt_level14_cumulative)
y_bacterias_pt = (bacterias_pt_level1_cumulative, bacterias_pt_level2_cumulative, bacterias_pt_level3_cumulative, bacterias_pt_level4_cumulative, bacterias_pt_level7_cumulative, bacterias_pt_level12_cumulative, bacterias_pt_level14_cumulative)
y_archaeas_pt = (archaeas_pt_level1_cumulative, archaeas_pt_level2_cumulative, archaeas_pt_level3_cumulative, archaeas_pt_level4_cumulative, archaeas_pt_level7_cumulative, archaeas_pt_level12_cumulative, archaeas_pt_level14_cumulative)
y_fungi_pt = (fungi_pt_level1_cumulative, fungi_pt_level2_cumulative, fungi_pt_level3_cumulative, fungi_pt_level4_cumulative, fungi_pt_level7_cumulative, fungi_pt_level12_cumulative, fungi_pt_level14_cumulative)
y_plant_pt = (plant_pt_level1_cumulative, plant_pt_level2_cumulative, plant_pt_level3_cumulative, plant_pt_level4_cumulative, plant_pt_level7_cumulative, plant_pt_level12_cumulative, plant_pt_level14_cumulative)
y_protozoa_pt = (protozoa_pt_level1_cumulative, protozoa_pt_level2_cumulative, protozoa_pt_level3_cumulative, protozoa_pt_level4_cumulative, protozoa_pt_level7_cumulative, protozoa_pt_level12_cumulative, protozoa_pt_level14_cumulative)
y_mito_pt = (mito_pt_level1_cumulative, mito_pt_level2_cumulative, mito_pt_level3_cumulative, mito_pt_level4_cumulative, mito_pt_level7_cumulative, mito_pt_level12_cumulative, mito_pt_level14_cumulative)
y_plastid_pt = (plastid_pt_level1_cumulative, plastid_pt_level2_cumulative, plastid_pt_level3_cumulative, plastid_pt_level4_cumulative, plastid_pt_level7_cumulative, plastid_pt_level12_cumulative, plastid_pt_level14_cumulative)


# plt.figure(figsize=(6,2), dpi=80)

# plt.bar(X,y_virus_pt, align='center', width=0.8)
# #plt.bar(x,y_virus_nm, align='center', width=0.4)
# plt.xlabel('Levels')
# plt.ylabel('Cumulative Bytes')
# plt.title("Level study")

# axes = plt.gca()
# #axes.set_ylim([98000000,100000000])
# axes.set_ylim([25000000,26000000])

# #for i in range(len(y)):
# #    plt.hlines(y[i],0,x[i])
# plt.show()


# plt.figure()
# plt.title("Virus Cumulative Benchmark")
# plt.ylabel('Cumulative Bytes')
# plt.xlabel('Levels')

'''
    Draw the plots for the different domains
'''

## Virus Level Analysis
fig_vir= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_virus_nm, align='center', color="blue")
ax1.set_ylim([98000000,100000000])
ax1.title.set_text('Virus Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_virus_pt, align='center', color="orange")
ax2.set_ylim([25000000,26000000])
ax2.title.set_text('Virus Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Virus.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Bacterias Level Analysis

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Archaea Level Analysis
fig_arch= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_archaeas_nm, align='center', color="blue")
ax1.set_ylim([730000000,740000000])
ax1.title.set_text('Archaea Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_archaeas_pt, align='center', color="orange")
ax2.set_ylim([210000000,220000000])
ax2.title.set_text('Archaea Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Archaea.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Fungi Level Analysis
fig_fungi= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_fungi_nm, align='center', color="blue")
ax1.set_ylim([260000000,280000000])
ax1.title.set_text('Fungi Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_fungi_pt, align='center', color="orange")
ax2.set_ylim([78000000,83000000])
ax2.title.set_text('Fungi Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Fungi.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Plant Level Analysis
fig_plant= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_plant_nm, align='center', color="blue")
ax1.set_ylim([120000000,140000000])
ax1.title.set_text('Plant Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_plant_pt, align='center', color="orange")
ax2.set_ylim([45000000,50000000])
ax2.title.set_text('Plant Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Plant.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Protozoa Level Analysis
fig_prot= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_protozoa_nm, align='center', color="blue")
ax1.set_ylim([49000000,57000000])
ax1.title.set_text('Protozoa Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_protozoa_pt, align='center', color="orange")
ax2.set_ylim([20000000,24000000])
ax2.title.set_text('Protozoa Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Protozoa.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Mitochondrial Level Analysis
fig_mito= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_mito_nm, align='center', color="blue")
ax1.set_ylim([69000000,71000000])
ax1.title.set_text('Mitochondrial Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_mito_pt, align='center', color="orange")
ax2.set_ylim([10000000,12000000])
ax2.title.set_text('Mitochondrial Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Mitochondrial.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Plastid Level Analysis
fig_mito= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_plastid_nm, align='center', color="blue")
ax1.set_ylim([175000000,180000000])
ax1.title.set_text('Plastid Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_plastid_pt, align='center', color="orange")
ax2.set_ylim([30000000,30500000])
ax2.title.set_text('Plastid Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Plastid.png')