'''
    Usage: python3 src/find_greatest_length.py DNA
           python3 src/find_greatest_length.py AA
'''

import sys, csv, math
import os, os.path
import numpy as np

Type = str(sys.argv[1]) # DNA | AA

maxvaluefound_length_DNA = 0
maxvaluefound_length_AA = 0

## Virus Values
virus_nm_length, virus_pt_length = [], []

## Bacterias Values
bacteria_nm_length, bacteria_pt_length = [], []

## Archaeas Values
archaea_nm_length, archaea_pt_length = [], []

## Fungis Values
fungi_nm_length, fungi_pt_length = [], []

## Plants Values
plant_nm_length, plant_pt_length = [], []

## Protozoas Values
protozoa_nm_length, protozoa_pt_length = [], []

## Mitochondrials Values
mito_nm_length, mito_pt_length = [], []

## Platids Values
plastid_nm_length, plastid_pt_length = [], []

def virusData(virus_nm_length, virus_pt_length):
    # Viral CSV files
    with open('Analysis/LengthSeq/Virus/length_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            virus_nm_length.append(row[0].split("\t")[1])
            virus_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    virus_nm_length = np.array(virus_nm_length[1:], dtype=float)
    virus_pt_length = np.array(virus_pt_length[1:], dtype=float)

    return virus_nm_length, virus_pt_length

def bacteriaData(bacteria_nm_length, bacteria_pt_length):
    # Bacterias CSV files
    with open('Analysis/LengthSeq/Bacteria/length_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            bacteria_nm_length.append(row[0].split("\t")[1])
            bacteria_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    bacteria_nm_length = np.array(bacteria_nm_length[1:], dtype=float)
    bacteria_pt_length = np.array(bacteria_pt_length[1:], dtype=float)

    return bacteria_nm_length, bacteria_pt_length

def archaeaData(archaea_nm_length, archaea_pt_length):
    # Archaeas CSV files
    with open('Analysis/LengthSeq/Archaea/length_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            archaea_nm_length.append(row[0].split("\t")[1])
            archaea_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    archaea_nm_length = np.array(archaea_nm_length[1:], dtype=float)
    archaea_pt_length = np.array(archaea_pt_length[1:], dtype=float)

    return archaea_nm_length, archaea_pt_length

def fungiData(fungi_nm_length, fungi_pt_length):
    # Fungis CSV files
    with open('Analysis/LengthSeq/Fungi/length_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            fungi_nm_length.append(row[0].split("\t")[1])
            fungi_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    fungi_nm_length = np.array(fungi_nm_length[1:], dtype=float)
    fungi_pt_length = np.array(fungi_pt_length[1:], dtype=float)

    return fungi_nm_length, fungi_pt_length

def plantData(plant_nm_length, plant_pt_length):
    # Plants CSV files
    with open('Analysis/LengthSeq/Plant/length_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plant_nm_length.append(row[0].split("\t")[1])
            plant_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    plant_nm_length = np.array(plant_nm_length[1:], dtype=float)
    plant_pt_length = np.array(plant_pt_length[1:], dtype=float)

    return plant_nm_length, plant_pt_length

def protozoaData(protozoa_nm_length, protozoa_pt_length):
    # Protozoa CSV files
    with open('Analysis/LengthSeq/Protozoa/length_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            protozoa_nm_length.append(row[0].split("\t")[1])
            protozoa_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    protozoa_nm_length = np.array(protozoa_nm_length[1:], dtype=float)
    protozoa_pt_length = np.array(protozoa_pt_length[1:], dtype=float)

    return protozoa_nm_length, protozoa_pt_length

def mitoData(mito_nm_length, mito_pt_length):
    # Mitochondrial CSV files
    with open('Analysis/LengthSeq/Mitochondrial/length_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            mito_nm_length.append(row[0].split("\t")[1])
            mito_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    mito_nm_length = np.array(mito_nm_length[1:], dtype=float)
    mito_pt_length = np.array(mito_pt_length[1:], dtype=float)

    return mito_nm_length, mito_pt_length

def plastidData(plastid_nm_length, plastid_pt_length):
    # Plastid CSV files
    with open('Analysis/LengthSeq/Plastid/length_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            plastid_nm_length.append(row[0].split("\t")[1])
            plastid_pt_length.append(row[0].split("\t")[2])

    # Data conversion and header removal
    plastid_nm_length = np.array(plastid_nm_length[1:], dtype=float)
    plastid_pt_length = np.array(plastid_pt_length[1:], dtype=float)

    return plastid_nm_length, plastid_pt_length

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

virus_nm_length, virus_pt_length = virusData(virus_nm_length, virus_pt_length)
bacteria_nm_length, bacteria_pt_length = bacteriaData(bacteria_nm_length, bacteria_pt_length)
archaea_nm_length, archaea_pt_length = archaeaData(archaea_nm_length, archaea_pt_length)
fungi_nm_length, fungi_pt_length = fungiData(fungi_nm_length, fungi_pt_length)
plant_nm_length, plant_pt_length = plantData(plant_nm_length, plant_pt_length)
protozoa_nm_length, protozoa_pt_length = protozoaData(protozoa_nm_length, protozoa_pt_length)
mito_nm_length, mito_pt_length = mitoData(mito_nm_length, mito_pt_length)
plastid_nm_length, plastid_pt_length = plastidData(plastid_nm_length, plastid_pt_length)

virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length, maxvaluefound_length_DNA = normalizeLengths(virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length)
virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length, maxvaluefound_length_AA = normalizeLengths(virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length)

if Type == "DNA":
    print(int(maxvaluefound_length_DNA))
elif Type == "AA":
    print(int(maxvaluefound_length_AA))