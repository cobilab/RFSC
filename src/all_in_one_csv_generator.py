import csv
import os, os.path
import numpy as np
import sys

L3 = 8                          # Level 3 (position 8 of the csv)

samples = int(sys.argv[1])       # Number of samples used in the csv                       

domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]

maxvaluefound_length_DNA = 0
maxvaluefound_length_AA = 0

dir_path = os.path.dirname(os.path.realpath(__file__))  # Current directory path

## Virus Values
virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length = [], [], [], [], []

## Bacterias Values
bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length = [], [], [], [], []

## Archaeas Values
archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length = [], [], [], [], []

## Fungis Values
fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length = [], [], [], [], []

## Plants Values
plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length = [], [], [], [], []

## Protozoas Values
protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length = [], [], [], [], []

## Mitochondrials Values
mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length = [], [], [], [], []

## Platids Values
plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length = [], [], [], [], []

def virusData(virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Virus/length_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                virus_nm_length.append(row[0].split("\t")[1])
                virus_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    virus_nm_result3 = np.array(virus_nm_result3[1:], dtype=float)
    virus_pt_result3 = np.array(virus_pt_result3[1:], dtype=float)
    virus_gc_content = np.array(virus_gc_content[1:], dtype=float)
    virus_nm_length = np.array(virus_nm_length[1:], dtype=float)
    virus_pt_length = np.array(virus_pt_length[1:], dtype=float)

    return virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length

def bacteriaData(bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length):
    
    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Bacteria/length_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                bacteria_nm_length.append(row[0].split("\t")[1])
                bacteria_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    bacteria_nm_result3 = np.array(bacteria_nm_result3[1:], dtype=float)
    bacteria_pt_result3 = np.array(bacteria_pt_result3[1:], dtype=float)
    bacteria_gc_content = np.array(bacteria_gc_content[1:], dtype=float)
    bacteria_nm_length = np.array(bacteria_nm_length[1:], dtype=float)
    bacteria_pt_length = np.array(bacteria_pt_length[1:], dtype=float)

    return bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length

def archaeaData(archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Archaea/length_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                archaea_nm_length.append(row[0].split("\t")[1])
                archaea_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    archaea_nm_result3 = np.array(archaea_nm_result3[1:], dtype=float)
    archaea_pt_result3 = np.array(archaea_pt_result3[1:], dtype=float)
    archaea_gc_content = np.array(archaea_gc_content[1:], dtype=float)
    archaea_nm_length = np.array(archaea_nm_length[1:], dtype=float)
    archaea_pt_length = np.array(archaea_pt_length[1:], dtype=float)

    return archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length

def fungiData(fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Fungi/length_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                fungi_nm_length.append(row[0].split("\t")[1])
                fungi_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    fungi_nm_result3 = np.array(fungi_nm_result3[1:], dtype=float)
    fungi_pt_result3 = np.array(fungi_pt_result3[1:], dtype=float)
    fungi_gc_content = np.array(fungi_gc_content[1:], dtype=float)
    fungi_nm_length = np.array(fungi_nm_length[1:], dtype=float)
    fungi_pt_length = np.array(fungi_pt_length[1:], dtype=float)

    return fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length

def plantData(plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Plant/length_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                plant_nm_length.append(row[0].split("\t")[1])
                plant_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    plant_nm_result3 = np.array(plant_nm_result3[1:], dtype=float)
    plant_pt_result3 = np.array(plant_pt_result3[1:], dtype=float)
    plant_gc_content = np.array(plant_gc_content[1:], dtype=float)
    plant_nm_length = np.array(plant_nm_length[1:], dtype=float)
    plant_pt_length = np.array(plant_pt_length[1:], dtype=float)

    return plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length

def protozoaData(protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                protozoa_nm_length.append(row[0].split("\t")[1])
                protozoa_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    protozoa_nm_result3 = np.array(protozoa_nm_result3[1:], dtype=float)
    protozoa_pt_result3 = np.array(protozoa_pt_result3[1:], dtype=float)
    protozoa_gc_content = np.array(protozoa_gc_content[1:], dtype=float)
    protozoa_nm_length = np.array(protozoa_nm_length[1:], dtype=float)
    protozoa_pt_length = np.array(protozoa_pt_length[1:], dtype=float)

    return protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length

def mitoData(mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Mitochondrial/length_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                mito_nm_length.append(row[0].split("\t")[1])
                mito_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    mito_nm_result3 = np.array(mito_nm_result3[1:], dtype=float)
    mito_pt_result3 = np.array(mito_pt_result3[1:], dtype=float)
    mito_gc_content = np.array(mito_gc_content[1:], dtype=float)
    mito_nm_length = np.array(mito_nm_length[1:], dtype=float)
    mito_pt_length = np.array(mito_pt_length[1:], dtype=float)

    return mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length

def plastidData(plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

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

    with open('Analysis/LengthSeq/Plastid/length_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if countLen <= samples:
                plastid_nm_length.append(row[0].split("\t")[1])
                plastid_pt_length.append(row[0].split("\t")[2])
                countLen = countLen + 1

    # Data conversion and header removal
    plastid_nm_result3 = np.array(plastid_nm_result3[1:], dtype=float)
    plastid_pt_result3 = np.array(plastid_pt_result3[1:], dtype=float)
    plastid_gc_content = np.array(plastid_gc_content[1:], dtype=float)
    plastid_nm_length = np.array(plastid_nm_length[1:], dtype=float)
    plastid_pt_length = np.array(plastid_pt_length[1:], dtype=float)

    return plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length

def normalizeLengths(virus_length, bacteria_length, archaea_length, fungi_lenth, plant_length, protozoa_length, mito_length, plastid_length):
    
    maxvalues = []

    maxvalues.append(max(virus_length))
    maxvalues.append(max(bacteria_length))
    maxvalues.append(max(archaea_length))
    maxvalues.append(max(fungi_lenth))
    maxvalues.append(max(plant_length))
    maxvalues.append(max(protozoa_length))
    maxvalues.append(max(mito_length))
    maxvalues.append(max(plastid_length))

    maxvaluefound = max(maxvalues)

    virus_length[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in virus_length]
    bacteria_length[:] = [float("{:.10f}".format(x/maxvaluefound))for x in bacteria_length]
    archaea_length[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in archaea_length]
    fungi_lenth[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in fungi_lenth]
    plant_length[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in plant_length]
    protozoa_length[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in protozoa_length]
    mito_length[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in mito_length]
    plastid_length[:] = [float("{:.10f}".format(x/maxvaluefound)) for x in plastid_length]

    return virus_length, bacteria_length, archaea_length, fungi_lenth, plant_length, protozoa_length, mito_length, plastid_length, maxvaluefound

virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length = virusData(virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length)
bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length = bacteriaData(bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length)
archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length = archaeaData(archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length)
fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length = fungiData(fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length)
plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length = plantData(plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length)
protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length = protozoaData(protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length)
mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length = mitoData(mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length)
plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length = plastidData(plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length)

virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length, maxvaluefound_length_DNA = normalizeLengths(virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length)
virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length, maxvaluefound_length_AA = normalizeLengths(virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length)

with open('Analysis/KNN/Domains.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Domain", "DNA", "AA", "GC", "L_DNA", "L_AA"])
    #for domain in domains:
    for x in range(samples):
        if x == len(virus_nm_result3):
            break
        writer.writerow([domains[0], virus_nm_result3[x], virus_pt_result3[x], virus_gc_content[x], float("{:.16f}".format(virus_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(virus_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(bacteria_nm_result3):
            break
        writer.writerow([domains[1], bacteria_nm_result3[x], bacteria_pt_result3[x], bacteria_gc_content[x], float("{:.16f}".format(bacteria_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(bacteria_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(archaea_nm_result3):
            break
        writer.writerow([domains[2], archaea_nm_result3[x], archaea_pt_result3[x], archaea_gc_content[x], float("{:.16f}".format(archaea_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(archaea_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(fungi_nm_result3):
            break
        writer.writerow([domains[3], fungi_nm_result3[x], fungi_pt_result3[x], fungi_gc_content[x], float("{:.16f}".format(fungi_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(fungi_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(plant_nm_result3):
            break
        writer.writerow([domains[4], plant_nm_result3[x], plant_pt_result3[x], plant_gc_content[x], float("{:.16f}".format(plant_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(plant_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(protozoa_nm_result3):
            break
        writer.writerow([domains[5], protozoa_nm_result3[x], protozoa_pt_result3[x], protozoa_gc_content[x], float("{:.16f}".format(protozoa_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(protozoa_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(mito_nm_result3):
            break
        writer.writerow([domains[6], mito_nm_result3[x], mito_pt_result3[x], mito_gc_content[x], float("{:.16f}".format(mito_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(mito_pt_length[x]/maxvaluefound_length_AA))])
    for x in range(samples):
        if x == len(plastid_nm_result3):
            break
        writer.writerow([domains[7], plastid_nm_result3[x], plastid_pt_result3[x], plastid_gc_content[x], float("{:.16f}".format(plastid_nm_length[x]/maxvaluefound_length_DNA)), float("{:.16f}".format(plastid_pt_length[x]/maxvaluefound_length_AA))])