import csv
import os, os.path
import numpy as np

L3 = 8              # Level 3 (position 8 of the csv)

samples = 500       # Number of samples used in the csv                       

domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]

dir_path = os.path.dirname(os.path.realpath(__file__))  # Current directory path

## Virus Values
virus_nm_result3, virus_pt_result3, virus_gc_content = [], [], []

## Bacterias Values
bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content = [], [], []

## Archaeas Values
archaea_nm_result3, archaea_pt_result3, archaea_gc_content = [], [], []

## Fungis Values
fungi_nm_result3, fungi_pt_result3, fungi_gc_content = [], [], []

## Plants Values
plant_nm_result3, plant_pt_result3, plant_gc_content = [], [], []

## Protozoas Values
protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content = [], [], []

## Mitochondrials Values
mito_nm_result3, mito_pt_result3, mito_gc_content = [], [], []

## Platids Values
plastid_nm_result3, plastid_pt_result3, plastid_gc_content = [], [], []

def virusData(virus_nm_result3, virus_pt_result3, virus_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Viral CSV files
    with open('Analysis/GeCo/Virus/geco3_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                virus_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Virus/ac2_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                virus_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Virus/gc_content_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                virus_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    virus_nm_result3 = np.array(virus_nm_result3[1:], dtype=float)
    virus_pt_result3 = np.array(virus_pt_result3[1:], dtype=float)
    virus_gc_content = np.array(virus_gc_content[1:], dtype=float)

    return virus_nm_result3, virus_pt_result3, virus_gc_content

def bacteriaData(bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content):
    
    countNM, countPT, countGC = 0, 0, 0

    # Bacterias CSV files
    with open('Analysis/GeCo/Bacteria/geco3_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                bacteria_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Bacteria/ac2_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                bacteria_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Bacteria/gc_content_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                bacteria_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    bacteria_nm_result3 = np.array(bacteria_nm_result3[1:], dtype=float)
    bacteria_pt_result3 = np.array(bacteria_pt_result3[1:], dtype=float)
    bacteria_gc_content = np.array(bacteria_gc_content[1:], dtype=float)

    return bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content

def archaeaData(archaea_nm_result3, archaea_pt_result3, archaea_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Archaeas CSV files
    with open('Analysis/GeCo/Archaea/geco3_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                archaea_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Archaea/ac2_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                archaea_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Archaea/gc_content_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                archaea_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    archaea_nm_result3 = np.array(archaea_nm_result3[1:], dtype=float)
    archaea_pt_result3 = np.array(archaea_pt_result3[1:], dtype=float)
    archaea_gc_content = np.array(archaea_gc_content[1:], dtype=float)

    return archaea_nm_result3, archaea_pt_result3, archaea_gc_content

def fungiData(fungi_nm_result3, fungi_pt_result3, fungi_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Fungis CSV files
    with open('Analysis/GeCo/Fungi/geco3_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                fungi_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Fungi/ac2_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                fungi_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Fungi/gc_content_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                fungi_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    fungi_nm_result3 = np.array(fungi_nm_result3[1:], dtype=float)
    fungi_pt_result3 = np.array(fungi_pt_result3[1:], dtype=float)
    fungi_gc_content = np.array(fungi_gc_content[1:], dtype=float)

    return fungi_nm_result3, fungi_pt_result3, fungi_gc_content

def plantData(plant_nm_result3, plant_pt_result3, plant_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Plants CSV files
    with open('Analysis/GeCo/Plant/geco3_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                plant_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Plant/ac2_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                plant_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Plant/gc_content_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                plant_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    plant_nm_result3 = np.array(plant_nm_result3[1:], dtype=float)
    plant_pt_result3 = np.array(plant_pt_result3[1:], dtype=float)
    plant_gc_content = np.array(plant_gc_content[1:], dtype=float)

    return plant_nm_result3, plant_pt_result3, plant_gc_content

def protozoaData(protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Protozoa CSV files
    with open('Analysis/GeCo/Protozoa/geco3_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                protozoa_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Protozoa/ac2_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                protozoa_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                protozoa_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    protozoa_nm_result3 = np.array(protozoa_nm_result3[1:], dtype=float)
    protozoa_pt_result3 = np.array(protozoa_pt_result3[1:], dtype=float)
    protozoa_gc_content = np.array(protozoa_gc_content[1:], dtype=float)

    return protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content

def mitoData(mito_nm_result3, mito_pt_result3, mito_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Mitochondrial CSV files
    with open('Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                mito_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Mitochondrial/ac2_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                mito_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Mitochondrial/gc_content_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                mito_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    mito_nm_result3 = np.array(mito_nm_result3[1:], dtype=float)
    mito_pt_result3 = np.array(mito_pt_result3[1:], dtype=float)
    mito_gc_content = np.array(mito_gc_content[1:], dtype=float)

    return mito_nm_result3, mito_pt_result3, mito_gc_content

def plastidData(plastid_nm_result3, plastid_pt_result3, plastid_gc_content):

    countNM, countPT, countGC = 0, 0, 0

    # Plastid CSV files
    with open('Analysis/GeCo/Plastid/geco3_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countNM <= samples:
                plastid_nm_result3.append("0"+row[0].split("\t")[L3])
                countNM = countNM + 1

    with open('Analysis/AC/Plastid/ac2_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countPT <= samples:
                plastid_pt_result3.append(row[0].split("\t")[L3])
                countPT = countPT + 1

    with open('Analysis/GCcontent/Plastid/gc_content_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countGC <= samples:
                plastid_gc_content.append("0"+row[0].split("\t")[3])
                countGC = countGC + 1

    # Data conversion and header removal
    plastid_nm_result3 = np.array(plastid_nm_result3[1:], dtype=float)
    plastid_pt_result3 = np.array(plastid_pt_result3[1:], dtype=float)
    plastid_gc_content = np.array(plastid_gc_content[1:], dtype=float)

    return plastid_nm_result3, plastid_pt_result3, plastid_gc_content


virus_nm_result3, virus_pt_result3, virus_gc_content = virusData(virus_nm_result3, virus_pt_result3, virus_gc_content)
bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content = bacteriaData(bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content)
archaea_nm_result3, archaea_pt_result3, archaea_gc_content = archaeaData(archaea_nm_result3, archaea_pt_result3, archaea_gc_content)
fungi_nm_result3, fungi_pt_result3, fungi_gc_content = fungiData(fungi_nm_result3, fungi_pt_result3, fungi_gc_content)
plant_nm_result3, plant_pt_result3, plant_gc_content = plantData(plant_nm_result3, plant_pt_result3, plant_gc_content)
protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content = protozoaData(protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content)
mito_nm_result3, mito_pt_result3, mito_gc_content = mitoData(mito_nm_result3, mito_pt_result3, mito_gc_content)
plastid_nm_result3, plastid_pt_result3, plastid_gc_content = plastidData(plastid_nm_result3, plastid_pt_result3, plastid_gc_content)

with open('Analysis/SVM/Domains.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Domain", "DNA", "AA", "GC"])
    for domain in domains:
        for x in range(samples):
            writer.writerow([domain, virus_nm_result3[x], virus_pt_result3[x], virus_gc_content[x]])