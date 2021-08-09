import csv
import os, os.path
import numpy as np
import random
import sys

L3 = 8                                  # Level 3 (position 8 of the csv)                    

train_partition = float(sys.argv[1])    # 0.75

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

    # Viral CSV files
    with open('Analysis/GeCo/Virus/geco3_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            virus_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Virus/ac2_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            virus_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Virus/gc_content_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            virus_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Virus/length_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            virus_nm_length.append(row[0].split("\t")[1])
            virus_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    virus_nm_result3 = np.array(virus_nm_result3[1:], dtype=float)
    virus_pt_result3 = np.array(virus_pt_result3[1:], dtype=float)
    virus_gc_content = np.array(virus_gc_content[1:], dtype=float)
    virus_nm_length = np.array(virus_nm_length[1:], dtype=float)
    virus_pt_length = np.array(virus_pt_length[1:], dtype=float)

    training_idx_viral = []
    for i in range(0,int(len(virus_nm_result3) * train_partition)):
        n = random.randint(1,len(virus_nm_result3))
        while n in training_idx_viral:
            n = random.randint(1,len(virus_nm_result3))
        training_idx_viral.append(n)

    return virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length, training_idx_viral

def bacteriaData(bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length):
    
    # Bacterias CSV files
    with open('Analysis/GeCo/Bacteria/geco3_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            bacteria_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Bacteria/ac2_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            bacteria_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Bacteria/gc_content_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            bacteria_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Bacteria/length_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            bacteria_nm_length.append(row[0].split("\t")[1])
            bacteria_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    bacteria_nm_result3 = np.array(bacteria_nm_result3[1:], dtype=float)
    bacteria_pt_result3 = np.array(bacteria_pt_result3[1:], dtype=float)
    bacteria_gc_content = np.array(bacteria_gc_content[1:], dtype=float)
    bacteria_nm_length = np.array(bacteria_nm_length[1:], dtype=float)
    bacteria_pt_length = np.array(bacteria_pt_length[1:], dtype=float)

    training_idx_bacteria = []
    for i in range(0,int(len(bacteria_nm_result3) * train_partition)):
        n = random.randint(1,len(bacteria_nm_result3))
        while n in training_idx_bacteria:
            n = random.randint(1,len(bacteria_nm_result3))
        training_idx_bacteria.append(n)

    return bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length, training_idx_bacteria

def archaeaData(archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length):

    # Archaeas CSV files
    with open('Analysis/GeCo/Archaea/geco3_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            archaea_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Archaea/ac2_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            archaea_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Archaea/gc_content_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            archaea_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Archaea/length_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            archaea_nm_length.append(row[0].split("\t")[1])
            archaea_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    archaea_nm_result3 = np.array(archaea_nm_result3[1:], dtype=float)
    archaea_pt_result3 = np.array(archaea_pt_result3[1:], dtype=float)
    archaea_gc_content = np.array(archaea_gc_content[1:], dtype=float)
    archaea_nm_length = np.array(archaea_nm_length[1:], dtype=float)
    archaea_pt_length = np.array(archaea_pt_length[1:], dtype=float)

    training_idx_archaea = []
    for i in range(0,int(len(archaea_nm_result3) * train_partition)):
        n = random.randint(1,len(archaea_nm_result3))
        while n in training_idx_archaea:
            n = random.randint(1,len(archaea_nm_result3))
        training_idx_archaea.append(n)

    return archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length, training_idx_archaea

def fungiData(fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length):

    # Fungis CSV files
    with open('Analysis/GeCo/Fungi/geco3_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            fungi_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Fungi/ac2_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            fungi_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Fungi/gc_content_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            fungi_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Fungi/length_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            fungi_nm_length.append(row[0].split("\t")[1])
            fungi_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    fungi_nm_result3 = np.array(fungi_nm_result3[1:], dtype=float)
    fungi_pt_result3 = np.array(fungi_pt_result3[1:], dtype=float)
    fungi_gc_content = np.array(fungi_gc_content[1:], dtype=float)
    fungi_nm_length = np.array(fungi_nm_length[1:], dtype=float)
    fungi_pt_length = np.array(fungi_pt_length[1:], dtype=float)

    training_idx_fungi = []
    for i in range(0,int(len(fungi_nm_result3) * train_partition)):
        n = random.randint(1,len(fungi_nm_result3))
        while n in training_idx_fungi:
            n = random.randint(1,len(fungi_nm_result3))
        training_idx_fungi.append(n)

    return fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length, training_idx_fungi

def plantData(plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length):

    # Plants CSV files
    with open('Analysis/GeCo/Plant/geco3_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plant_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Plant/ac2_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plant_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Plant/gc_content_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plant_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Plant/length_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plant_nm_length.append(row[0].split("\t")[1])
            plant_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    plant_nm_result3 = np.array(plant_nm_result3[1:], dtype=float)
    plant_pt_result3 = np.array(plant_pt_result3[1:], dtype=float)
    plant_gc_content = np.array(plant_gc_content[1:], dtype=float)
    plant_nm_length = np.array(plant_nm_length[1:], dtype=float)
    plant_pt_length = np.array(plant_pt_length[1:], dtype=float)

    training_idx_plant = []
    for i in range(0,int(len(plant_nm_result3) * train_partition)):
        n = random.randint(1,len(plant_nm_result3))
        while n in training_idx_plant:
            n = random.randint(1,len(plant_nm_result3))
        training_idx_plant.append(n)

    return plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length, training_idx_plant

def protozoaData(protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length):

    # Protozoa CSV files
    with open('Analysis/GeCo/Protozoa/geco3_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            protozoa_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Protozoa/ac2_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            protozoa_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            protozoa_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            protozoa_nm_length.append(row[0].split("\t")[1])
            protozoa_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    protozoa_nm_result3 = np.array(protozoa_nm_result3[1:], dtype=float)
    protozoa_pt_result3 = np.array(protozoa_pt_result3[1:], dtype=float)
    protozoa_gc_content = np.array(protozoa_gc_content[1:], dtype=float)
    protozoa_nm_length = np.array(protozoa_nm_length[1:], dtype=float)
    protozoa_pt_length = np.array(protozoa_pt_length[1:], dtype=float)

    training_idx_protozoa = []
    for i in range(0,int(len(protozoa_nm_result3) * train_partition)):
        n = random.randint(1,len(protozoa_nm_result3))
        while n in training_idx_protozoa:
            n = random.randint(1,len(protozoa_nm_result3))
        training_idx_protozoa.append(n)

    return protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length, training_idx_protozoa

def mitoData(mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length):

    # Mitochondrial CSV files
    with open('Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            mito_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Mitochondrial/ac2_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            mito_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Mitochondrial/gc_content_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            mito_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Mitochondrial/length_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            mito_nm_length.append(row[0].split("\t")[1])
            mito_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    mito_nm_result3 = np.array(mito_nm_result3[1:], dtype=float)
    mito_pt_result3 = np.array(mito_pt_result3[1:], dtype=float)
    mito_gc_content = np.array(mito_gc_content[1:], dtype=float)
    mito_nm_length = np.array(mito_nm_length[1:], dtype=float)
    mito_pt_length = np.array(mito_pt_length[1:], dtype=float)

    training_idx_mito = []
    for i in range(0,int(len(mito_nm_result3) * train_partition)):
        n = random.randint(1,len(mito_nm_result3))
        while n in training_idx_mito:
            n = random.randint(1,len(mito_nm_result3))
        training_idx_mito.append(n)

    return mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length, training_idx_mito

def plastidData(plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length):

    # Plastid CSV files
    with open('Analysis/GeCo/Plastid/geco3_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plastid_nm_result3.append("0"+row[0].split("\t")[L3])

    with open('Analysis/AC/Plastid/ac2_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plastid_pt_result3.append(row[0].split("\t")[L3])

    with open('Analysis/GCcontent/Plastid/gc_content_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plastid_gc_content.append("0"+row[0].split("\t")[3])

    with open('Analysis/LengthSeq/Plastid/length_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plastid_nm_length.append(row[0].split("\t")[1])
            plastid_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    plastid_nm_result3 = np.array(plastid_nm_result3[1:], dtype=float)
    plastid_pt_result3 = np.array(plastid_pt_result3[1:], dtype=float)
    plastid_gc_content = np.array(plastid_gc_content[1:], dtype=float)
    plastid_nm_length = np.array(plastid_nm_length[1:], dtype=float)
    plastid_pt_length = np.array(plastid_pt_length[1:], dtype=float)

    training_idx_plastid = []
    for i in range(0,int(len(plastid_nm_result3) * train_partition)):
        n = random.randint(1,len(plastid_nm_result3))
        while n in training_idx_plastid:
            n = random.randint(1,len(plastid_nm_result3))
        training_idx_plastid.append(n)

    return plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length, training_idx_plastid

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

virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length, training_idx_viral = virusData(virus_nm_result3, virus_pt_result3, virus_gc_content, virus_nm_length, virus_pt_length)
bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length, training_idx_bacteria = bacteriaData(bacteria_nm_result3, bacteria_pt_result3, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length)
archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length, training_idx_archaea = archaeaData(archaea_nm_result3, archaea_pt_result3, archaea_gc_content, archaea_nm_length, archaea_pt_length)
fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length, training_idx_fungi = fungiData(fungi_nm_result3, fungi_pt_result3, fungi_gc_content, fungi_nm_length, fungi_pt_length)
plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length, training_idx_plant = plantData(plant_nm_result3, plant_pt_result3, plant_gc_content, plant_nm_length, plant_pt_length)
protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length, training_idx_protozoa = protozoaData(protozoa_nm_result3, protozoa_pt_result3, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length)
mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length, training_idx_mito = mitoData(mito_nm_result3, mito_pt_result3, mito_gc_content, mito_nm_length, mito_pt_length)
plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length, training_idx_plastid = plastidData(plastid_nm_result3, plastid_pt_result3, plastid_gc_content, plastid_nm_length, plastid_pt_length)

virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length, maxvaluefound_length_DNA = normalizeLengths(virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length)
virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length, maxvaluefound_length_AA = normalizeLengths(virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length)

with open('Analysis/KNN/Train.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Domain", "DNA", "AA", "GC", "L_DNA", "L_AA"])
    #for domain in domains:
    for x in range(0,len(virus_nm_result3)):
        if x in training_idx_viral:
            writer.writerow([domains[0], virus_nm_result3[x], virus_pt_result3[x], virus_gc_content[x], virus_nm_length[x], virus_pt_length[x]])
    for x in range(0,len(bacteria_nm_result3)):
        if x in training_idx_bacteria:
            writer.writerow([domains[1], bacteria_nm_result3[x], bacteria_pt_result3[x], bacteria_gc_content[x], bacteria_nm_length[x], bacteria_pt_length[x]])
    for x in range(0,len(archaea_nm_result3)):
        if x in training_idx_archaea:
            writer.writerow([domains[2], archaea_nm_result3[x], archaea_pt_result3[x], archaea_gc_content[x], archaea_nm_length[x], archaea_pt_length[x]])
    for x in range(0,len(fungi_nm_result3)):
        if x in training_idx_fungi:
            writer.writerow([domains[3], fungi_nm_result3[x], fungi_pt_result3[x], fungi_gc_content[x], fungi_nm_length[x], fungi_pt_length[x]])
    for x in range(0,len(plant_nm_result3)):
        if x in training_idx_plant:
            writer.writerow([domains[4], plant_nm_result3[x], plant_pt_result3[x], plant_gc_content[x], plant_nm_length[x], plant_pt_length[x]])
    for x in range(0, len(protozoa_nm_result3)):
        if x in training_idx_protozoa:
            writer.writerow([domains[5], protozoa_nm_result3[x], protozoa_pt_result3[x], protozoa_gc_content[x], protozoa_nm_length[x], protozoa_pt_length[x]])
    for x in range(0, len(mito_nm_result3)):
        if x in training_idx_mito:
            writer.writerow([domains[6], mito_nm_result3[x], mito_pt_result3[x], mito_gc_content[x], mito_nm_length[x], mito_pt_length[x]])
    for x in range(0, len(plastid_nm_result3)):
        if x in training_idx_plastid:
            writer.writerow([domains[7], plastid_nm_result3[x], plastid_pt_result3[x], plastid_gc_content[x], plastid_nm_length[x], plastid_pt_length[x]])

with open('Analysis/KNN/Test.csv', 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(["Domain", "DNA", "AA", "GC", "L_DNA", "L_AA"])
    #for domain in domains:
    for x in range(0,len(virus_nm_result3)):
        if x not in training_idx_viral:
            writer.writerow([domains[0], virus_nm_result3[x], virus_pt_result3[x], virus_gc_content[x], virus_nm_length[x], virus_pt_length[x]])
    for x in range(0,len(bacteria_nm_result3)):
        if x not in training_idx_bacteria:
            writer.writerow([domains[1], bacteria_nm_result3[x], bacteria_pt_result3[x], bacteria_gc_content[x], bacteria_nm_length[x], bacteria_pt_length[x]])
    for x in range(0,len(archaea_nm_result3)):
        if x not in training_idx_archaea:
            writer.writerow([domains[2], archaea_nm_result3[x], archaea_pt_result3[x], archaea_gc_content[x], archaea_nm_length[x], archaea_pt_length[x]])
    for x in range(0,len(fungi_nm_result3)):
        if x not in training_idx_fungi:
            writer.writerow([domains[3], fungi_nm_result3[x], fungi_pt_result3[x], fungi_gc_content[x], fungi_nm_length[x], fungi_pt_length[x]])
    for x in range(0,len(plant_nm_result3)):
        if x not in training_idx_plant:
            writer.writerow([domains[4], plant_nm_result3[x], plant_pt_result3[x], plant_gc_content[x], plant_nm_length[x], plant_pt_length[x]])
    for x in range(0, len(protozoa_nm_result3)):
        if x not in training_idx_protozoa:
            writer.writerow([domains[5], protozoa_nm_result3[x], protozoa_pt_result3[x], protozoa_gc_content[x], protozoa_nm_length[x], protozoa_pt_length[x]])
    for x in range(0, len(mito_nm_result3)):
        if x not in training_idx_mito:
            writer.writerow([domains[6], mito_nm_result3[x], mito_pt_result3[x], mito_gc_content[x], mito_nm_length[x], mito_pt_length[x]])
    for x in range(0, len(plastid_nm_result3)):
        if x not in training_idx_plastid:
            writer.writerow([domains[7], plastid_nm_result3[x], plastid_pt_result3[x], plastid_gc_content[x], plastid_nm_length[x], plastid_pt_length[x]])