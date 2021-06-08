"""
    Generation of a Benchmark for Level os Compression study

    x: Levels of study
    y: Cumulative Sizes of the sequences of the given domain 

    Run: python3 levelsBenchmark.py
"""

from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv

## Virus Values
virus_nm_size1, virus_nm_size2, virus_nm_size3, virus_nm_size4, virus_nm_size7 = [], [], [], [], []
virus_pt_size1, virus_pt_size2, virus_pt_size3, virus_pt_size4, virus_pt_size7 = [], [], [], [], []

## Bacterias Values
bacterias_nm_size1, bacterias_nm_size2, bacterias_nm_size3, bacterias_nm_size4, bacterias_nm_size7 = [], [], [], [], []
bacterias_pt_size1, bacterias_pt_size2, bacterias_pt_size3, bacterias_pt_size4, bacterias_pt_size7 = [], [], [], [], []

## Archaeas Values
archaeas_nm_size1, archaeas_nm_size2, archaeas_nm_size3, archaeas_nm_size4, archaeas_nm_size7, = [], [], [], [], []
archaeas_pt_size1, archaeas_pt_size2, archaeas_pt_size3, archaeas_pt_size4, archaeas_pt_size7, = [], [], [], [], []

## Fungi Values
fungi_nm_size1, fungi_nm_size2, fungi_nm_size3, fungi_nm_size4, fungi_nm_size7 = [], [], [], [], []
fungi_pt_size1, fungi_pt_size2, fungi_pt_size3, fungi_pt_size4, fungi_pt_size7 = [], [], [], [], []

## Plant Values
plant_nm_size1, plant_nm_size2, plant_nm_size3, plant_nm_size4, plant_nm_size7 = [], [], [], [], []
plant_pt_size1, plant_pt_size2, plant_pt_size3, plant_pt_size4, plant_pt_size7 = [], [], [], [], []

## Protozoa Values
protozoa_nm_size1, protozoa_nm_size2, protozoa_nm_size3, protozoa_nm_size4, protozoa_nm_size7 = [], [], [], [], []
protozoa_pt_size1, protozoa_pt_size2, protozoa_pt_size3, protozoa_pt_size4, protozoa_pt_size7 = [], [], [], [], []

## Mitochondrial Values
mito_nm_size1, mito_nm_size2, mito_nm_size3, mito_nm_size4, mito_nm_size7 = [], [], [], [], []
mito_pt_size1, mito_pt_size2, mito_pt_size3, mito_pt_size4, mito_pt_size7 = [], [], [], [], []

## Plastid Values
plastid_nm_size1, plastid_nm_size2, plastid_nm_size3, plastid_nm_size4, plastid_nm_size7 = [], [], [], [], []
plastid_pt_size1, plastid_pt_size2, plastid_pt_size3, plastid_pt_size4, plastid_pt_size7 = [], [], [], [], []

# Virus CSV files
with open('Analysis/GeCo/Virus/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_nm_size1.append(row[0].split("\t")[1])
        virus_nm_size2.append(row[0].split("\t")[2])
        virus_nm_size3.append(row[0].split("\t")[3])
        virus_nm_size4.append(row[0].split("\t")[4])
        virus_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Virus/ac2_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_pt_size1.append(row[0].split("\t")[1])
        virus_pt_size2.append(row[0].split("\t")[2])
        virus_pt_size3.append(row[0].split("\t")[3])
        virus_pt_size4.append(row[0].split("\t")[4])
        virus_pt_size7.append(row[0].split("\t")[5])

# Bacterias CSV files
with open('Analysis/GeCo/Bacteria/geco3_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacterias_nm_size1.append(row[0].split("\t")[1])
        bacterias_nm_size2.append(row[0].split("\t")[2])
        bacterias_nm_size3.append(row[0].split("\t")[3])
        bacterias_nm_size4.append(row[0].split("\t")[4])
        bacterias_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Bacteria/ac2_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacterias_pt_size1.append(row[0].split("\t")[1])
        bacterias_pt_size2.append(row[0].split("\t")[2])
        bacterias_pt_size3.append(row[0].split("\t")[3])
        bacterias_pt_size4.append(row[0].split("\t")[4])
        bacterias_pt_size7.append(row[0].split("\t")[5])


# Archaeas CSV files
with open('Analysis/GeCo/Archaea/geco3_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaeas_nm_size1.append(row[0].split("\t")[1])
        archaeas_nm_size2.append(row[0].split("\t")[2])
        archaeas_nm_size3.append(row[0].split("\t")[3])
        archaeas_nm_size4.append(row[0].split("\t")[4])
        archaeas_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Archaea/ac2_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaeas_pt_size1.append(row[0].split("\t")[1])
        archaeas_pt_size2.append(row[0].split("\t")[2])
        archaeas_pt_size3.append(row[0].split("\t")[3])
        archaeas_pt_size4.append(row[0].split("\t")[4])
        archaeas_pt_size7.append(row[0].split("\t")[5])

# Fungis CSV files
with open('Analysis/GeCo/Fungi/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_nm_size1.append(row[0].split("\t")[1])
        fungi_nm_size2.append(row[0].split("\t")[2])
        fungi_nm_size3.append(row[0].split("\t")[3])
        fungi_nm_size4.append(row[0].split("\t")[4])
        fungi_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Fungi/ac2_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_pt_size1.append(row[0].split("\t")[1])
        fungi_pt_size2.append(row[0].split("\t")[2])
        fungi_pt_size3.append(row[0].split("\t")[3])
        fungi_pt_size4.append(row[0].split("\t")[4])
        fungi_pt_size7.append(row[0].split("\t")[5])

# Plants CSV files
with open('Analysis/GeCo/Plant/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_nm_size1.append(row[0].split("\t")[1])
        plant_nm_size2.append(row[0].split("\t")[2])
        plant_nm_size3.append(row[0].split("\t")[3])
        plant_nm_size4.append(row[0].split("\t")[4])
        plant_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Plant/ac2_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_pt_size1.append(row[0].split("\t")[1])
        plant_pt_size2.append(row[0].split("\t")[2])
        plant_pt_size3.append(row[0].split("\t")[3])
        plant_pt_size4.append(row[0].split("\t")[4])
        plant_pt_size7.append(row[0].split("\t")[5])

# Protozoa CSV files
with open('Analysis/GeCo/Protozoa/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_nm_size1.append(row[0].split("\t")[1])
        protozoa_nm_size2.append(row[0].split("\t")[2])
        protozoa_nm_size3.append(row[0].split("\t")[3])
        protozoa_nm_size4.append(row[0].split("\t")[4])
        protozoa_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Protozoa/ac2_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_pt_size1.append(row[0].split("\t")[1])
        protozoa_pt_size2.append(row[0].split("\t")[2])
        protozoa_pt_size3.append(row[0].split("\t")[3])
        protozoa_pt_size4.append(row[0].split("\t")[4])
        protozoa_pt_size7.append(row[0].split("\t")[5])

# Mitochondrial CSV files
with open('Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_nm_size1.append(row[0].split("\t")[1])
        mito_nm_size2.append(row[0].split("\t")[2])
        mito_nm_size3.append(row[0].split("\t")[3])
        mito_nm_size4.append(row[0].split("\t")[4])
        mito_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Mitochondrial/ac2_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_pt_size1.append(row[0].split("\t")[1])
        mito_pt_size2.append(row[0].split("\t")[2])
        mito_pt_size3.append(row[0].split("\t")[3])
        mito_pt_size4.append(row[0].split("\t")[4])
        mito_pt_size7.append(row[0].split("\t")[5])

# Plastid CSV files
with open('Analysis/GeCo/Plastid/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_nm_size1.append(row[0].split("\t")[1])
        plastid_nm_size2.append(row[0].split("\t")[2])
        plastid_nm_size3.append(row[0].split("\t")[3])
        plastid_nm_size4.append(row[0].split("\t")[4])
        plastid_nm_size7.append(row[0].split("\t")[5])

with open('Analysis/AC/Plastid/ac2_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_pt_size1.append(row[0].split("\t")[1])
        plastid_pt_size2.append(row[0].split("\t")[2])
        plastid_pt_size3.append(row[0].split("\t")[3])
        plastid_pt_size4.append(row[0].split("\t")[4])
        plastid_pt_size7.append(row[0].split("\t")[5])

'''
    Convert sizes from string to int
'''

## Virus conversion
virus_nm_size1 = np.array(virus_nm_size1[1:], dtype=int)
virus_nm_size2 = np.array(virus_nm_size2[1:], dtype=int)
virus_nm_size3 = np.array(virus_nm_size3[1:], dtype=int)
virus_nm_size4 = np.array(virus_nm_size4[1:], dtype=int)
virus_nm_size7 = np.array(virus_nm_size7[1:], dtype=int)

virus_pt_size1 = np.array(virus_pt_size1[1:], dtype=int)
virus_pt_size2 = np.array(virus_pt_size2[1:], dtype=int)
virus_pt_size3 = np.array(virus_pt_size3[1:], dtype=int)
virus_pt_size4 = np.array(virus_pt_size4[1:], dtype=int)
virus_pt_size7 = np.array(virus_pt_size7[1:], dtype=int)

## Bacterias conversion
bacterias_nm_size1 = np.array(bacterias_nm_size1[1:], dtype=int)
bacterias_nm_size2 = np.array(bacterias_nm_size2[1:], dtype=int)
bacterias_nm_size3 = np.array(bacterias_nm_size3[1:], dtype=int)
bacterias_nm_size4 = np.array(bacterias_nm_size4[1:], dtype=int)
bacterias_nm_size7 = np.array(bacterias_nm_size7[1:], dtype=int)

bacterias_pt_size1 = np.array(bacterias_pt_size1[1:], dtype=int)
bacterias_pt_size2 = np.array(bacterias_pt_size2[1:], dtype=int)
bacterias_pt_size3 = np.array(bacterias_pt_size3[1:], dtype=int)
bacterias_pt_size4 = np.array(bacterias_pt_size4[1:], dtype=int)
bacterias_pt_size7 = np.array(bacterias_pt_size7[1:], dtype=int)

## Archaeas conversion
archaeas_nm_size1 = np.array(archaeas_nm_size1[1:], dtype=int)
archaeas_nm_size2 = np.array(archaeas_nm_size2[1:], dtype=int)
archaeas_nm_size3 = np.array(archaeas_nm_size3[1:], dtype=int)
archaeas_nm_size4 = np.array(archaeas_nm_size4[1:], dtype=int)
archaeas_nm_size7 = np.array(archaeas_nm_size7[1:], dtype=int)

archaeas_pt_size1 = np.array(archaeas_pt_size1[1:], dtype=int)
archaeas_pt_size2 = np.array(archaeas_pt_size2[1:], dtype=int)
archaeas_pt_size3 = np.array(archaeas_pt_size3[1:], dtype=int)
archaeas_pt_size4 = np.array(archaeas_pt_size4[1:], dtype=int)
archaeas_pt_size7 = np.array(archaeas_pt_size7[1:], dtype=int)

## Fungis conversion
fungi_nm_size1 = np.array(fungi_nm_size1[1:], dtype=int)
fungi_nm_size2 = np.array(fungi_nm_size2[1:], dtype=int)
fungi_nm_size3 = np.array(fungi_nm_size3[1:], dtype=int)
fungi_nm_size4 = np.array(fungi_nm_size4[1:], dtype=int)
fungi_nm_size7 = np.array(fungi_nm_size7[1:], dtype=int)

fungi_pt_size1 = np.array(fungi_pt_size1[1:], dtype=int)
fungi_pt_size2 = np.array(fungi_pt_size2[1:], dtype=int)
fungi_pt_size3 = np.array(fungi_pt_size3[1:], dtype=int)
fungi_pt_size4 = np.array(fungi_pt_size4[1:], dtype=int)
fungi_pt_size7 = np.array(fungi_pt_size7[1:], dtype=int)

## Plant conversion
plant_nm_size1 = np.array(plant_nm_size1[1:], dtype=int)
plant_nm_size2 = np.array(plant_nm_size2[1:], dtype=int)
plant_nm_size3 = np.array(plant_nm_size3[1:], dtype=int)
plant_nm_size4 = np.array(plant_nm_size4[1:], dtype=int)
plant_nm_size7 = np.array(plant_nm_size7[1:], dtype=int)

plant_pt_size1 = np.array(plant_pt_size1[1:], dtype=int)
plant_pt_size2 = np.array(plant_pt_size2[1:], dtype=int)
plant_pt_size3 = np.array(plant_pt_size3[1:], dtype=int)
plant_pt_size4 = np.array(plant_pt_size4[1:], dtype=int)
plant_pt_size7 = np.array(plant_pt_size7[1:], dtype=int)

## Protozoa conversion
protozoa_nm_size1 = np.array(protozoa_nm_size1[1:], dtype=int)
protozoa_nm_size2 = np.array(protozoa_nm_size2[1:], dtype=int)
protozoa_nm_size3 = np.array(protozoa_nm_size3[1:], dtype=int)
protozoa_nm_size4 = np.array(protozoa_nm_size4[1:], dtype=int)
protozoa_nm_size7 = np.array(protozoa_nm_size7[1:], dtype=int)

protozoa_pt_size1 = np.array(protozoa_pt_size1[1:], dtype=int)
protozoa_pt_size2 = np.array(protozoa_pt_size2[1:], dtype=int)
protozoa_pt_size3 = np.array(protozoa_pt_size3[1:], dtype=int)
protozoa_pt_size4 = np.array(protozoa_pt_size4[1:], dtype=int)
protozoa_pt_size7 = np.array(protozoa_pt_size7[1:], dtype=int)

## Mitochondrial conversion
mito_nm_size1 = np.array(mito_nm_size1[1:], dtype=int)
mito_nm_size2 = np.array(mito_nm_size2[1:], dtype=int)
mito_nm_size3 = np.array(mito_nm_size3[1:], dtype=int)
mito_nm_size4 = np.array(mito_nm_size4[1:], dtype=int)
mito_nm_size7 = np.array(mito_nm_size7[1:], dtype=int)

mito_pt_size1 = np.array(mito_pt_size1[1:], dtype=int)
mito_pt_size2 = np.array(mito_pt_size2[1:], dtype=int)
mito_pt_size3 = np.array(mito_pt_size3[1:], dtype=int)
mito_pt_size4 = np.array(mito_pt_size4[1:], dtype=int)
mito_pt_size7 = np.array(mito_pt_size7[1:], dtype=int)

## Plastid conversion
plastid_nm_size1 = np.array(plastid_nm_size1[1:], dtype=int)
plastid_nm_size2 = np.array(plastid_nm_size2[1:], dtype=int)
plastid_nm_size3 = np.array(plastid_nm_size3[1:], dtype=int)
plastid_nm_size4 = np.array(plastid_nm_size4[1:], dtype=int)
plastid_nm_size7 = np.array(plastid_nm_size7[1:], dtype=int)

plastid_pt_size1 = np.array(plastid_pt_size1[1:], dtype=int)
plastid_pt_size2 = np.array(plastid_pt_size2[1:], dtype=int)
plastid_pt_size3 = np.array(plastid_pt_size3[1:], dtype=int)
plastid_pt_size4 = np.array(plastid_pt_size4[1:], dtype=int)
plastid_pt_size7 = np.array(plastid_pt_size7[1:], dtype=int)

'''
    Calculate cumulative sizes for every level of each domain
'''

# sum(Virus)
virus_nm_level1_cumulative = sum(virus_nm_size1)
virus_nm_level2_cumulative = sum(virus_nm_size2)
virus_nm_level3_cumulative = sum(virus_nm_size3)
virus_nm_level4_cumulative = sum(virus_nm_size4)
virus_nm_level7_cumulative = sum(virus_nm_size7)

virus_pt_level1_cumulative = sum(virus_pt_size1)
virus_pt_level2_cumulative = sum(virus_pt_size2)
virus_pt_level3_cumulative = sum(virus_pt_size3)
virus_pt_level4_cumulative = sum(virus_pt_size4)
virus_pt_level7_cumulative = sum(virus_pt_size7)

# sum(Bacteria)
bacterias_nm_level1_cumulative = sum(bacterias_nm_size1)
bacterias_nm_level2_cumulative = sum(bacterias_nm_size2)
bacterias_nm_level3_cumulative = sum(bacterias_nm_size3)
bacterias_nm_level4_cumulative = sum(bacterias_nm_size4)
bacterias_nm_level7_cumulative = sum(bacterias_nm_size7)

bacterias_pt_level1_cumulative = sum(bacterias_pt_size1)
bacterias_pt_level2_cumulative = sum(bacterias_pt_size2)
bacterias_pt_level3_cumulative = sum(bacterias_pt_size3)
bacterias_pt_level4_cumulative = sum(bacterias_pt_size4)
bacterias_pt_level7_cumulative = sum(bacterias_pt_size7)

# sum(Archaea)
archaeas_nm_level1_cumulative = sum(archaeas_nm_size1)
archaeas_nm_level2_cumulative = sum(archaeas_nm_size2)
archaeas_nm_level3_cumulative = sum(archaeas_nm_size3)
archaeas_nm_level4_cumulative = sum(archaeas_nm_size4)
archaeas_nm_level7_cumulative = sum(archaeas_nm_size7)

archaeas_pt_level1_cumulative = sum(archaeas_pt_size1)
archaeas_pt_level2_cumulative = sum(archaeas_pt_size2)
archaeas_pt_level3_cumulative = sum(archaeas_pt_size3)
archaeas_pt_level4_cumulative = sum(archaeas_pt_size4)
archaeas_pt_level7_cumulative = sum(archaeas_pt_size7)

# sum(Fungi)
fungi_nm_level1_cumulative = sum(fungi_nm_size1)
fungi_nm_level2_cumulative = sum(fungi_nm_size2)
fungi_nm_level3_cumulative = sum(fungi_nm_size3)
fungi_nm_level4_cumulative = sum(fungi_nm_size4)
fungi_nm_level7_cumulative = sum(fungi_nm_size7)

fungi_pt_level1_cumulative = sum(fungi_pt_size1)
fungi_pt_level2_cumulative = sum(fungi_pt_size2)
fungi_pt_level3_cumulative = sum(fungi_pt_size3)
fungi_pt_level4_cumulative = sum(fungi_pt_size4)
fungi_pt_level7_cumulative = sum(fungi_pt_size7)

# sum(Plant)
plant_nm_level1_cumulative = sum(plant_nm_size1)
plant_nm_level2_cumulative = sum(plant_nm_size2)
plant_nm_level3_cumulative = sum(plant_nm_size3)
plant_nm_level4_cumulative = sum(plant_nm_size4)
plant_nm_level7_cumulative = sum(plant_nm_size7)

plant_pt_level1_cumulative = sum(plant_pt_size1)
plant_pt_level2_cumulative = sum(plant_pt_size2)
plant_pt_level3_cumulative = sum(plant_pt_size3)
plant_pt_level4_cumulative = sum(plant_pt_size4)
plant_pt_level7_cumulative = sum(plant_pt_size7)

# sum(Protozoa)
protozoa_nm_level1_cumulative = sum(protozoa_nm_size1)
protozoa_nm_level2_cumulative = sum(protozoa_nm_size2)
protozoa_nm_level3_cumulative = sum(protozoa_nm_size3)
protozoa_nm_level4_cumulative = sum(protozoa_nm_size4)
protozoa_nm_level7_cumulative = sum(protozoa_nm_size7)

protozoa_pt_level1_cumulative = sum(protozoa_pt_size1)
protozoa_pt_level2_cumulative = sum(protozoa_pt_size2)
protozoa_pt_level3_cumulative = sum(protozoa_pt_size3)
protozoa_pt_level4_cumulative = sum(protozoa_pt_size4)
protozoa_pt_level7_cumulative = sum(protozoa_pt_size7)

# sum(Mitochondrial)
mito_nm_level1_cumulative = sum(mito_nm_size1)
mito_nm_level2_cumulative = sum(mito_nm_size2)
mito_nm_level3_cumulative = sum(mito_nm_size3)
mito_nm_level4_cumulative = sum(mito_nm_size4)
mito_nm_level7_cumulative = sum(mito_nm_size7)

mito_pt_level1_cumulative = sum(mito_pt_size1)
mito_pt_level2_cumulative = sum(mito_pt_size2)
mito_pt_level3_cumulative = sum(mito_pt_size3)
mito_pt_level4_cumulative = sum(mito_pt_size4)
mito_pt_level7_cumulative = sum(mito_pt_size7)

# sum(Plastid)
plastid_nm_level1_cumulative = sum(plastid_nm_size1)
plastid_nm_level2_cumulative = sum(plastid_nm_size2)
plastid_nm_level3_cumulative = sum(plastid_nm_size3)
plastid_nm_level4_cumulative = sum(plastid_nm_size4)
plastid_nm_level7_cumulative = sum(plastid_nm_size7)

plastid_pt_level1_cumulative = sum(plastid_pt_size1)
plastid_pt_level2_cumulative = sum(plastid_pt_size2)
plastid_pt_level3_cumulative = sum(plastid_pt_size3)
plastid_pt_level4_cumulative = sum(plastid_pt_size4)
plastid_pt_level7_cumulative = sum(plastid_pt_size7)


X = ['1', '2', '3', '4', '7']
y_virus_nm = (virus_nm_level1_cumulative, virus_nm_level2_cumulative, virus_nm_level3_cumulative, virus_nm_level4_cumulative, virus_nm_level7_cumulative)
y_bacterias_nm = (bacterias_nm_level1_cumulative, bacterias_nm_level2_cumulative, bacterias_nm_level3_cumulative, bacterias_nm_level4_cumulative, bacterias_nm_level7_cumulative)
y_archaeas_nm = (archaeas_nm_level1_cumulative, archaeas_nm_level2_cumulative, archaeas_nm_level3_cumulative, archaeas_nm_level4_cumulative, archaeas_nm_level7_cumulative)
y_fungi_nm = (fungi_nm_level1_cumulative, fungi_nm_level2_cumulative, fungi_nm_level3_cumulative, fungi_nm_level4_cumulative, fungi_nm_level7_cumulative)
y_plant_nm = (plant_nm_level1_cumulative, plant_nm_level2_cumulative, plant_nm_level3_cumulative, plant_nm_level4_cumulative, plant_nm_level7_cumulative)
y_protozoa_nm = (protozoa_nm_level1_cumulative, protozoa_nm_level2_cumulative, protozoa_nm_level3_cumulative, protozoa_nm_level4_cumulative, protozoa_nm_level7_cumulative)
y_mito_nm = (mito_nm_level1_cumulative, mito_nm_level2_cumulative, mito_nm_level3_cumulative, mito_nm_level4_cumulative, mito_nm_level7_cumulative)
y_plastid_nm = (plastid_nm_level1_cumulative, plastid_nm_level2_cumulative, plastid_nm_level3_cumulative, plastid_nm_level4_cumulative, plastid_nm_level7_cumulative)

y_virus_pt = (virus_pt_level1_cumulative, virus_pt_level2_cumulative, virus_pt_level3_cumulative, virus_pt_level4_cumulative, virus_pt_level7_cumulative)
y_bacterias_pt = (bacterias_pt_level1_cumulative, bacterias_pt_level2_cumulative, bacterias_pt_level3_cumulative, bacterias_pt_level4_cumulative, bacterias_pt_level7_cumulative)
y_archaeas_pt = (archaeas_pt_level1_cumulative, archaeas_pt_level2_cumulative, archaeas_pt_level3_cumulative, archaeas_pt_level4_cumulative, archaeas_pt_level7_cumulative)
y_fungi_pt = (fungi_pt_level1_cumulative, fungi_pt_level2_cumulative, fungi_pt_level3_cumulative, fungi_pt_level4_cumulative, fungi_pt_level7_cumulative)
y_plant_pt = (plant_pt_level1_cumulative, plant_pt_level2_cumulative, plant_pt_level3_cumulative, plant_pt_level4_cumulative, plant_pt_level7_cumulative)
y_protozoa_pt = (protozoa_pt_level1_cumulative, protozoa_pt_level2_cumulative, protozoa_pt_level3_cumulative, protozoa_pt_level4_cumulative, protozoa_pt_level7_cumulative)
y_mito_pt = (mito_pt_level1_cumulative, mito_pt_level2_cumulative, mito_pt_level3_cumulative, mito_pt_level4_cumulative, mito_pt_level7_cumulative)
y_plastid_pt = (plastid_pt_level1_cumulative, plastid_pt_level2_cumulative, plastid_pt_level3_cumulative, plastid_pt_level4_cumulative, plastid_pt_level7_cumulative)


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
ax2.set_ylim([240000000,250000000])
ax2.title.set_text('Virus Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Virus.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Bacterias Level Analysis
fig_bac= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_bacterias_nm, align='center', color="blue")
ax1.set_ylim([1300000000,1400000000]) # Change
ax1.title.set_text('Bacteria Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_bacterias_pt, align='center', color="orange")
ax2.set_ylim([1900000000,1950000000]) # Change
ax2.title.set_text('Bacteria Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Bacteria.png')

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

plt.clf() # Clear all previous plots

## Archaea Level Analysis
fig_arch= plt.subplots(1,2,figsize=(6,4))
ax1 = plt.subplot(221)
ax2 = plt.subplot(222)


ax1.bar(X, y_archaeas_nm, align='center', color="blue")
ax1.set_ylim([750000000,785000000])
ax1.title.set_text('Archaea Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_archaeas_pt, align='center', color="orange")
ax2.set_ylim([2200000000,2250000000])
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
ax1.set_ylim([285000000,295000000])
ax1.title.set_text('Fungi Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_fungi_pt, align='center', color="orange")
ax2.set_ylim([900000000,920000000])
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
ax1.set_ylim([140000000,155000000])
ax1.title.set_text('Plant Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_plant_pt, align='center', color="orange")
ax2.set_ylim([545000000,620000000])
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
ax1.set_ylim([62000000,63500000])
ax1.title.set_text('Protozoa Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_protozoa_pt, align='center', color="orange")
ax2.set_ylim([260000000,280000000])
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
ax1.set_ylim([71000000,72000000])
ax1.title.set_text('Mitochondrial Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_mito_pt, align='center', color="orange")
ax2.set_ylim([120000000,135000000])
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
ax1.set_ylim([180000000,190000000])
ax1.title.set_text('Plastid Nucleotide')
ax1.set_ylabel('Cumulative Bytes')
ax1.set_xlabel('Levels')

ax2.bar(X, y_plastid_pt, align='center', color="orange")
ax2.set_ylim([360000000,370000000])
ax2.title.set_text('Plastid Protein')
ax2.set_xlabel('Levels')

plt.savefig('Analysis/Level_Analysis/Plastid.png')