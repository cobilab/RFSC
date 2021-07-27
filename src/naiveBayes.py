import sys, csv, math
import os, os.path
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from statistics import NormalDist

L1 = 6                              # Level 1 (position 6 of the csv)
L3 = 8                              # Level 3 (position 8 of the csv)
L7 = 10                             # Level 7 (position 10 of the csv)

N = int(sys.argv[1])                # Number of types (Domains)

trainDatabase = float(sys.argv[2])  # Percentage of the (current) domain database used for training the model (80%)

crossValidation_Min_Lim = float(sys.argv[3])    # Inferior limit of the block used for testing
crossValidation_Max_Lim = float(sys.argv[4])    # Superior limit of the block used for testing

sample_DNA_1 = float(sys.argv[5])       # Normalize Compression Rate of the DNA for the Input Sequence (Level 1)
sample_AA_1 = float(sys.argv[6])        # Normalize Compression Rate of the AA for the Input Sequence (Level 1)
sample_DNA_3 = float(sys.argv[7])       # Normalize Compression Rate of the DNA for the Input Sequence (Level 3)
sample_AA_3 = float(sys.argv[8])        # Normalize Compression Rate of the AA for the Input Sequence (Level 3)
sample_DNA_7 = float(sys.argv[9])       # Normalize Compression Rate of the DNA for the Input Sequence (Level 7)
sample_AA_7 = float(sys.argv[10])       # Normalize Compression Rate of the AA for the Input Sequence (Level 7)
gc_percent = float(sys.argv[11])        # GC-Content Percentage ([0..1]) of the Input Sequence
sample_DNA_LEN = float(sys.argv[12])    # Length of the DNA
sample_AA_LEN = float(sys.argv[13])     # Length of the AA

preditors_required = str(sys.argv[14])

maxvaluefound_length_DNA = 0
maxvaluefound_length_AA = 0

dir_path = os.path.dirname(os.path.realpath(__file__))  # Current directory path

## Virus Values
virus_nm_result1, virus_pt_result1, virus_nm_result3, virus_pt_result3, virus_nm_result7, virus_pt_result7, virus_gc_content, virus_nm_length, virus_pt_length = [], [], [], [], [], [], [], [], []

## Bacterias Values
bacteria_nm_result1, bacteria_pt_result1, bacteria_nm_result3, bacteria_pt_result3, bacteria_nm_result7, bacteria_pt_result7, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length = [], [], [], [], [], [], [], [], []

## Archaeas Values
archaea_nm_result1, archaea_pt_result1, archaea_nm_result3, archaea_pt_result3, archaea_nm_result7, archaea_pt_result7, archaea_gc_content, archaea_nm_length, archaea_pt_length = [], [], [], [], [], [], [], [], []

## Fungis Values
fungi_nm_result1, fungi_pt_result1, fungi_nm_result3, fungi_pt_result3, fungi_nm_result7, fungi_pt_result7, fungi_gc_content, fungi_nm_length, fungi_pt_length = [], [], [], [], [], [], [], [], []

## Plants Values
plant_nm_result1, plant_pt_result1, plant_nm_result3, plant_pt_result3, plant_nm_result7, plant_pt_result7, plant_gc_content, plant_nm_length, plant_pt_length = [], [], [], [], [], [], [], [], []

## Protozoas Values
protozoa_nm_result1, protozoa_pt_result1, protozoa_nm_result3, protozoa_pt_result3, protozoa_nm_result7, protozoa_pt_result7, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length = [], [], [], [], [], [], [], [], []

## Mitochondrials Values
mito_nm_result1, mito_pt_result1, mito_nm_result3, mito_pt_result3, mito_nm_result7, mito_pt_result7, mito_gc_content, mito_nm_length, mito_pt_length = [], [], [], [], [], [], [], [], []

## Platids Values
plastid_nm_result1, plastid_pt_result1, plastid_nm_result3, plastid_pt_result3, plastid_nm_result7, plastid_pt_result7, plastid_gc_content, plastid_nm_length, plastid_pt_length = [], [], [], [], [], [], [], [], []

'''
    Generate Histogram Function
    Example: generateHist(virus_nm_result7, mu_virus_nm, std_virus_nm)
'''
def generateHist(data, mu, std):
    # Plot the histogram.
    plt.hist(data, bins=25, density=True, alpha=0.6, color='g')

    # Plot the PDF.
    xmin, xmax = plt.xlim()
    x = np.linspace(xmin, xmax, 100)
    p_virus_nm = norm.pdf(x, mu, std)

    plt.plot(x, p_virus_nm, 'k', linewidth=2)
    title = "Viral Nucleotides: \u03BC=%.2f,  \u03C3=%.2f" % (mu, std)
    plt.title(title)
    plt.ylabel("Frequency (%)")
    plt.xlabel("Normalized Dissimilarity Rate")

    plt.show()

def generateHist_GC():
    fig, ((ax0, ax1, ax2, ax3), (ax4, ax5, ax6, ax7)) = plt.subplots(nrows=2, ncols=4, figsize=(10,4))

    ax0.hist(virus_gc_content, bins=25, density=True, alpha=0.6, color='g')
    virus_xmin, virus_xmax = plt.xlim()
    x_virus = np.linspace(virus_xmin, virus_xmax, 100)
    p_virus_gc = norm.pdf(x_virus, mu_virus_gc, std_virus_gc)
    ax0.set_title("Viral GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_virus_gc, std_virus_gc))
    ax0.plot(x_virus, p_virus_gc, 'k', linewidth=2)

    ax1.hist(bacteria_gc_content, bins=25, density=True, alpha=0.6, color='g')
    bacteria_xmin, bacteria_xmax = plt.xlim()
    x_bacteria = np.linspace(bacteria_xmin, bacteria_xmax, 100)
    p_bacteria_gc = norm.pdf(x_bacteria, mu_bacteria_gc, std_bacteria_gc)
    ax1.set_title("Bacteria GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_bacteria_gc, std_bacteria_gc))
    ax1.plot(x_bacteria, p_bacteria_gc, 'k', linewidth=2)

    ax2.hist(archaea_gc_content, bins=25, density=True, alpha=0.6, color='g')
    archaea_xmin, archaea_xmax = plt.xlim()
    x_archaea = np.linspace(archaea_xmin, archaea_xmax, 100)
    p_archaea_gc = norm.pdf(x_archaea, mu_archaea_gc, std_archaea_gc)
    ax2.set_title("Archaea GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_archaea_gc, std_archaea_gc))
    ax2.plot(x_archaea, p_archaea_gc, 'k', linewidth=2)

    ax3.hist(fungi_gc_content, bins=25, density=True, alpha=0.6, color='g')
    fungi_xmin, fungi_xmax = plt.xlim()
    x_fungi = np.linspace(fungi_xmin, fungi_xmax, 100)
    p_fungi_gc = norm.pdf(x_fungi, mu_fungi_gc, std_fungi_gc)
    ax3.set_title("Fungi GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_fungi_gc, std_fungi_gc))
    ax3.plot(x_fungi, p_fungi_gc, 'k', linewidth=2)

    ax4.hist(plant_gc_content, bins=25, density=True, alpha=0.6, color='g')
    plant_xmin, plant_xmax = plt.xlim()
    x_plant = np.linspace(plant_xmin, plant_xmax, 100)
    p_plant_gc = norm.pdf(x_plant, mu_plant_gc, std_plant_gc)
    ax4.set_title("Plant GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_plant_gc, std_plant_gc))
    ax4.plot(x_plant, p_plant_gc, 'k', linewidth=2)

    ax5.hist(protozoa_gc_content, bins=25, density=True, alpha=0.6, color='g')
    protozoa_xmin, protozoa_xmax = plt.xlim()
    x_protozoa = np.linspace(protozoa_xmin, protozoa_xmax, 100)
    p_protozoa_gc = norm.pdf(x_protozoa, mu_protozoa_gc, std_protozoa_gc)
    ax5.set_title("Protozoa GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_protozoa_gc, std_protozoa_gc))
    ax5.plot(x_protozoa, p_protozoa_gc, 'k', linewidth=2)

    ax6.hist(mito_gc_content, bins=25, density=True, alpha=0.6, color='g')
    mito_xmin, mito_xmax = plt.xlim()
    x_mito = np.linspace(mito_xmin, mito_xmax, 100)
    p_mito_gc = norm.pdf(x_mito, mu_mito_gc, std_mito_gc)
    ax6.set_title("Protozoa GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_mito_gc, std_mito_gc))
    ax6.plot(x_mito, p_mito_gc, 'k', linewidth=2)

    ax7.hist(plastid_gc_content, bins=25, density=True, alpha=0.6, color='g')
    plastid_xmin, plastid_xmax = plt.xlim()
    x_plastid = np.linspace(plastid_xmin, plastid_xmax, 100)
    p_plastid_gc = norm.pdf(x_plastid, mu_plastid_gc, std_plastid_gc)
    ax7.set_title("Protozoa GC-Content: ")#\u03BC=%.2f,  \u03C3=%.2f" % (mu_plastid_gc, std_plastid_gc))
    ax7.plot(x_plastid, p_plastid_gc, 'k', linewidth=2)

    fig.tight_layout()
    plt.show()

def virusData(virus_nm_result1, virus_pt_result1, virus_nm_result3, virus_pt_result3, virus_nm_result7, virus_pt_result7, virus_gc_content, virus_nm_length, virus_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    virusFileNM = dir_path + '/../Analysis/GeCo/Virus/geco3_Viral.csv'
    virusFilePT = dir_path + '/../Analysis/AC/Virus/ac2_Viral.csv'
    virusFileGC = dir_path + '/../Analysis/GCcontent/Virus/gc_content_Viral.csv'
    virusFileLN = dir_path + '/../Analysis/LengthSeq/Virus/length_Viral.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(virusFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(virusFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(virusFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(virusFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(virusFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(virusFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(virusFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(virusFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(virusFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(virusFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(virusFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(virusFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Viral CSV files
    with open('Analysis/GeCo/Virus/geco3_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    virus_nm_result1.append("0"+row[0].split("\t")[L1])
                    virus_nm_result3.append("0"+row[0].split("\t")[L3])
                    virus_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    virus_nm_result1.append("0"+row[0].split("\t")[L1])
                    virus_nm_result3.append("0"+row[0].split("\t")[L3])
                    virus_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Virus/ac2_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    virus_pt_result1.append(row[0].split("\t")[L1])
                    virus_pt_result3.append(row[0].split("\t")[L3])
                    virus_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    virus_pt_result1.append(row[0].split("\t")[L1])
                    virus_pt_result3.append(row[0].split("\t")[L3])
                    virus_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Virus/gc_content_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    virus_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    virus_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Virus/length_Viral.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    virus_nm_length.append(row[0].split("\t")[1])
                    virus_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    virus_nm_length.append(row[0].split("\t")[1])
                    virus_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    virus_nm_result1 = np.array(virus_nm_result1[1:], dtype=float)
    virus_pt_result1 = np.array(virus_pt_result1[1:], dtype=float)
    virus_nm_result3 = np.array(virus_nm_result3[1:], dtype=float)
    virus_pt_result3 = np.array(virus_pt_result3[1:], dtype=float)
    virus_nm_result7 = np.array(virus_nm_result7[1:], dtype=float)
    virus_pt_result7 = np.array(virus_pt_result7[1:], dtype=float)
    virus_gc_content = np.array(virus_gc_content[1:], dtype=float)
    virus_nm_length = np.array(virus_nm_length[1:], dtype=float)
    virus_pt_length = np.array(virus_pt_length[1:], dtype=float)

    return virus_nm_result1, virus_pt_result1, virus_nm_result3, virus_pt_result3, virus_nm_result7, virus_pt_result7, virus_gc_content, virus_nm_length, virus_pt_length

def bacteriaData(bacteria_nm_result1, bacteria_pt_result1, bacteria_nm_result3, bacteria_pt_result3, bacteria_nm_result7, bacteria_pt_result7, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length):
    
    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    bacteriaFileNM = dir_path + '/../Analysis/GeCo/Bacteria/geco3_Bacteria.csv'
    bacteriaFilePT = dir_path + '/../Analysis/AC/Bacteria/ac2_Bacteria.csv'
    bacteriaFileGC = dir_path + '/../Analysis/GCcontent/Bacteria/gc_content_Bacteria.csv'
    bacteriaFileLN = dir_path + '/../Analysis/LengthSeq/Bacteria/length_Bacteria.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(bacteriaFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(bacteriaFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(bacteriaFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(bacteriaFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(bacteriaFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(bacteriaFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(bacteriaFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(bacteriaFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(bacteriaFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(bacteriaFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(bacteriaFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(bacteriaFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Bacterias CSV files
    with open('Analysis/GeCo/Bacteria/geco3_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    bacteria_nm_result1.append("0"+row[0].split("\t")[L1])
                    bacteria_nm_result3.append("0"+row[0].split("\t")[L3])
                    bacteria_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    bacteria_nm_result1.append("0"+row[0].split("\t")[L1])
                    bacteria_nm_result3.append("0"+row[0].split("\t")[L3])
                    bacteria_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Bacteria/ac2_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    bacteria_pt_result1.append(row[0].split("\t")[L1])
                    bacteria_pt_result3.append(row[0].split("\t")[L3])
                    bacteria_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    bacteria_pt_result1.append(row[0].split("\t")[L1])
                    bacteria_pt_result3.append(row[0].split("\t")[L3])
                    bacteria_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Bacteria/gc_content_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    bacteria_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    bacteria_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Bacteria/length_Bacteria.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    bacteria_nm_length.append(row[0].split("\t")[1])
                    bacteria_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    bacteria_nm_length.append(row[0].split("\t")[1])
                    bacteria_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    bacteria_nm_result1 = np.array(bacteria_nm_result1[1:], dtype=float)
    bacteria_pt_result1 = np.array(bacteria_pt_result1[1:], dtype=float)
    bacteria_nm_result3 = np.array(bacteria_nm_result3[1:], dtype=float)
    bacteria_pt_result3 = np.array(bacteria_pt_result3[1:], dtype=float)
    bacteria_nm_result7 = np.array(bacteria_nm_result7[1:], dtype=float)
    bacteria_pt_result7 = np.array(bacteria_pt_result7[1:], dtype=float)
    bacteria_gc_content = np.array(bacteria_gc_content[1:], dtype=float)
    bacteria_nm_length = np.array(bacteria_nm_length[1:], dtype=float)
    bacteria_pt_length = np.array(bacteria_pt_length[1:], dtype=float)

    return bacteria_nm_result1, bacteria_pt_result1, bacteria_nm_result3, bacteria_pt_result3, bacteria_nm_result7, bacteria_pt_result7, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length

def archaeaData(archaea_nm_result1, archaea_pt_result1, archaea_nm_result3, archaea_pt_result3, archaea_nm_result7, archaea_pt_result7, archaea_gc_content, archaea_nm_length, archaea_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    archaeaFileNM = dir_path + '/../Analysis/GeCo/Archaea/geco3_Archaea.csv'
    archaeaFilePT = dir_path + '/../Analysis/AC/Archaea/ac2_Archaea.csv'
    archaeaFileGC = dir_path + '/../Analysis/GCcontent/Archaea/gc_content_Archaea.csv'
    archaeaFileLN = dir_path + '/../Analysis/LengthSeq/Archaea/length_Archaea.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(archaeaFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(archaeaFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(archaeaFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(archaeaFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(archaeaFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(archaeaFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(archaeaFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(archaeaFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(archaeaFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(archaeaFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(archaeaFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(archaeaFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Archaeas CSV files
    with open('Analysis/GeCo/Archaea/geco3_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    archaea_nm_result1.append("0"+row[0].split("\t")[L1])
                    archaea_nm_result3.append("0"+row[0].split("\t")[L3])
                    archaea_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    archaea_nm_result1.append("0"+row[0].split("\t")[L1])
                    archaea_nm_result3.append("0"+row[0].split("\t")[L3])
                    archaea_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Archaea/ac2_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    archaea_pt_result1.append(row[0].split("\t")[L1])
                    archaea_pt_result3.append(row[0].split("\t")[L3])
                    archaea_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    archaea_pt_result1.append(row[0].split("\t")[L1])
                    archaea_pt_result3.append(row[0].split("\t")[L3])
                    archaea_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Archaea/gc_content_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    archaea_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    archaea_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Archaea/length_Archaea.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    archaea_nm_length.append(row[0].split("\t")[1])
                    archaea_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    archaea_nm_length.append(row[0].split("\t")[1])
                    archaea_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    archaea_nm_result1 = np.array(archaea_nm_result1[1:], dtype=float)
    archaea_pt_result1 = np.array(archaea_pt_result1[1:], dtype=float)
    archaea_nm_result3 = np.array(archaea_nm_result3[1:], dtype=float)
    archaea_pt_result3 = np.array(archaea_pt_result3[1:], dtype=float)
    archaea_nm_result7 = np.array(archaea_nm_result7[1:], dtype=float)
    archaea_pt_result7 = np.array(archaea_pt_result7[1:], dtype=float)
    archaea_gc_content = np.array(archaea_gc_content[1:], dtype=float)
    archaea_nm_length = np.array(archaea_nm_length[1:], dtype=float)
    archaea_pt_length = np.array(archaea_pt_length[1:], dtype=float)

    return archaea_nm_result1, archaea_pt_result1, archaea_nm_result3, archaea_pt_result3, archaea_nm_result7, archaea_pt_result7, archaea_gc_content, archaea_nm_length, archaea_pt_length

def fungiData(fungi_nm_result1, fungi_pt_result1, fungi_nm_result3, fungi_pt_result3, fungi_nm_result7, fungi_pt_result7, fungi_gc_content, fungi_nm_length, fungi_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    fungiFileNM = dir_path + '/../Analysis/GeCo/Fungi/geco3_Fungi.csv'
    fungiFilePT = dir_path + '/../Analysis/AC/Fungi/ac2_Fungi.csv'
    fungiFileGC = dir_path + '/../Analysis/GCcontent/Fungi/gc_content_Fungi.csv'
    fungiFileLN = dir_path + '/../Analysis/LengthSeq/Fungi/length_Fungi.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(fungiFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(fungiFileNM).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(fungiFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(fungiFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(fungiFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(fungiFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(fungiFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(fungiFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(fungiFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(fungiFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(fungiFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(fungiFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Fungis CSV files
    with open('Analysis/GeCo/Fungi/geco3_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    fungi_nm_result1.append("0"+row[0].split("\t")[L1])
                    fungi_nm_result3.append("0"+row[0].split("\t")[L3])
                    fungi_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    fungi_nm_result1.append("0"+row[0].split("\t")[L1])
                    fungi_nm_result3.append("0"+row[0].split("\t")[L3])
                    fungi_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Fungi/ac2_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    fungi_pt_result1.append(row[0].split("\t")[L1])
                    fungi_pt_result3.append(row[0].split("\t")[L3])
                    fungi_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    fungi_pt_result1.append(row[0].split("\t")[L1])
                    fungi_pt_result3.append(row[0].split("\t")[L3])
                    fungi_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Fungi/gc_content_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    fungi_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    fungi_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Fungi/length_Fungi.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    fungi_nm_length.append(row[0].split("\t")[1])
                    fungi_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    fungi_nm_length.append(row[0].split("\t")[1])
                    fungi_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    fungi_nm_result1 = np.array(fungi_nm_result1[1:], dtype=float)
    fungi_pt_result1 = np.array(fungi_pt_result1[1:], dtype=float)
    fungi_nm_result3 = np.array(fungi_nm_result3[1:], dtype=float)
    fungi_pt_result3 = np.array(fungi_pt_result3[1:], dtype=float)
    fungi_nm_result7 = np.array(fungi_nm_result7[1:], dtype=float)
    fungi_pt_result7 = np.array(fungi_pt_result7[1:], dtype=float)
    fungi_gc_content = np.array(fungi_gc_content[1:], dtype=float)
    fungi_nm_length = np.array(fungi_nm_length[1:], dtype=float)
    fungi_pt_length = np.array(fungi_pt_length[1:], dtype=float)

    return fungi_nm_result1, fungi_pt_result1, fungi_nm_result3, fungi_pt_result3, fungi_nm_result7, fungi_pt_result7, fungi_gc_content, fungi_nm_length, fungi_pt_length

def plantData(plant_nm_result1, plant_pt_result1, plant_nm_result3, plant_pt_result3, plant_nm_result7, plant_pt_result7, plant_gc_content, plant_nm_length, plant_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    plantFileNM = dir_path + '/../Analysis/GeCo/Plant/geco3_Plant.csv'
    plantFilePT = dir_path + '/../Analysis/AC/Plant/ac2_Plant.csv'
    plantFileGC = dir_path + '/../Analysis/GCcontent/Plant/gc_content_Plant.csv'
    plantFileLN = dir_path + '/../Analysis/LengthSeq/Plant/length_Plant.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(plantFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(plantFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(plantFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(plantFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(plantFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(plantFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(plantFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(plantFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(plantFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(plantFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(plantFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(plantFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Plants CSV files
    with open('Analysis/GeCo/Plant/geco3_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    plant_nm_result1.append("0"+row[0].split("\t")[L1])
                    plant_nm_result3.append("0"+row[0].split("\t")[L3])
                    plant_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    plant_nm_result1.append("0"+row[0].split("\t")[L1])
                    plant_nm_result3.append("0"+row[0].split("\t")[L3])
                    plant_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Plant/ac2_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    plant_pt_result1.append(row[0].split("\t")[L1])
                    plant_pt_result3.append(row[0].split("\t")[L3])
                    plant_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    plant_pt_result1.append(row[0].split("\t")[L1])
                    plant_pt_result3.append(row[0].split("\t")[L3])
                    plant_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Plant/gc_content_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    plant_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    plant_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Plant/length_Plant.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    plant_nm_length.append(row[0].split("\t")[1])
                    plant_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    plant_nm_length.append(row[0].split("\t")[1])
                    plant_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    plant_nm_result1 = np.array(plant_nm_result1[1:], dtype=float)
    plant_pt_result1 = np.array(plant_pt_result1[1:], dtype=float)
    plant_nm_result3 = np.array(plant_nm_result3[1:], dtype=float)
    plant_pt_result3 = np.array(plant_pt_result3[1:], dtype=float)
    plant_nm_result7 = np.array(plant_nm_result7[1:], dtype=float)
    plant_pt_result7 = np.array(plant_pt_result7[1:], dtype=float)
    plant_gc_content = np.array(plant_gc_content[1:], dtype=float)
    plant_nm_length = np.array(plant_nm_length[1:], dtype=float)
    plant_pt_length = np.array(plant_pt_length[1:], dtype=float)

    return plant_nm_result1, plant_pt_result1, plant_nm_result3, plant_pt_result3, plant_nm_result7, plant_pt_result7, plant_gc_content, plant_nm_length, plant_pt_length

def protozoaData(protozoa_nm_result1, protozoa_pt_result1, protozoa_nm_result3, protozoa_pt_result3, protozoa_nm_result7, protozoa_pt_result7, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    protozoaFileNM = dir_path + '/../Analysis/GeCo/Protozoa/geco3_Protozoa.csv'
    protozoaFilePT = dir_path + '/../Analysis/AC/Protozoa/ac2_Protozoa.csv'
    protozoaFileGC = dir_path + '/../Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv'
    protozoaFileLN = dir_path + '/../Analysis/LengthSeq/Protozoa/length_Protozoa.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(protozoaFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(protozoaFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(protozoaFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(protozoaFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(protozoaFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(protozoaFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(protozoaFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(protozoaFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(protozoaFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(protozoaFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(protozoaFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(protozoaFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Protozoa CSV files
    with open('Analysis/GeCo/Protozoa/geco3_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    protozoa_nm_result1.append("0"+row[0].split("\t")[L1])
                    protozoa_nm_result3.append("0"+row[0].split("\t")[L3])
                    protozoa_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    protozoa_nm_result1.append("0"+row[0].split("\t")[L1])
                    protozoa_nm_result3.append("0"+row[0].split("\t")[L3])
                    protozoa_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Protozoa/ac2_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    protozoa_pt_result1.append(row[0].split("\t")[L1])
                    protozoa_pt_result3.append(row[0].split("\t")[L3])
                    protozoa_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    protozoa_pt_result1.append(row[0].split("\t")[L1])
                    protozoa_pt_result3.append(row[0].split("\t")[L3])
                    protozoa_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    protozoa_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    protozoa_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Protozoa/length_Protozoa.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    protozoa_nm_length.append(row[0].split("\t")[1])
                    protozoa_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    protozoa_nm_length.append(row[0].split("\t")[1])
                    protozoa_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    protozoa_nm_result1 = np.array(protozoa_nm_result1[1:], dtype=float)
    protozoa_pt_result1 = np.array(protozoa_pt_result1[1:], dtype=float)
    protozoa_nm_result3 = np.array(protozoa_nm_result3[1:], dtype=float)
    protozoa_pt_result3 = np.array(protozoa_pt_result3[1:], dtype=float)
    protozoa_nm_result7 = np.array(protozoa_nm_result7[1:], dtype=float)
    protozoa_pt_result7 = np.array(protozoa_pt_result7[1:], dtype=float)
    protozoa_gc_content = np.array(protozoa_gc_content[1:], dtype=float)
    protozoa_nm_length = np.array(protozoa_nm_length[1:], dtype=float)
    protozoa_pt_length = np.array(protozoa_pt_length[1:], dtype=float)

    return protozoa_nm_result1, protozoa_pt_result1, protozoa_nm_result3, protozoa_pt_result3, protozoa_nm_result7, protozoa_pt_result7, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length

def mitoData(mito_nm_result1, mito_pt_result1, mito_nm_result3, mito_pt_result3, mito_nm_result7, mito_pt_result7, mito_gc_content, mito_nm_length, mito_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    mitoFileNM = dir_path + '/../Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv'
    mitoFilePT = dir_path + '/../Analysis/AC/Mitochondrial/ac2_Mitochondrial.csv'
    mitoFileGC = dir_path + '/../Analysis/GCcontent/Mitochondrial/gc_content_Mitochondrial.csv'
    mitoFileLN = dir_path + '/../Analysis/LengthSeq/Mitochondrial/length_Mitochondrial.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(mitoFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(mitoFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(mitoFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(mitoFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(mitoFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(mitoFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(mitoFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(mitoFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(mitoFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(mitoFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(mitoFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(mitoFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Mitochondrial CSV files
    with open('Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    mito_nm_result1.append("0"+row[0].split("\t")[L1])
                    mito_nm_result3.append("0"+row[0].split("\t")[L3])
                    mito_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    mito_nm_result1.append("0"+row[0].split("\t")[L1])
                    mito_nm_result3.append("0"+row[0].split("\t")[L3])
                    mito_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Mitochondrial/ac2_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    mito_pt_result1.append(row[0].split("\t")[L1])
                    mito_pt_result3.append(row[0].split("\t")[L3])
                    mito_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    mito_pt_result1.append(row[0].split("\t")[L1])
                    mito_pt_result3.append(row[0].split("\t")[L3])
                    mito_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Mitochondrial/gc_content_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    mito_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    mito_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Mitochondrial/length_Mitochondrial.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    mito_nm_length.append(row[0].split("\t")[1])
                    mito_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    mito_nm_length.append(row[0].split("\t")[1])
                    mito_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    mito_nm_result1 = np.array(mito_nm_result1[1:], dtype=float)
    mito_pt_result1 = np.array(mito_pt_result1[1:], dtype=float)
    mito_nm_result3 = np.array(mito_nm_result3[1:], dtype=float)
    mito_pt_result3 = np.array(mito_pt_result3[1:], dtype=float)
    mito_nm_result7 = np.array(mito_nm_result7[1:], dtype=float)
    mito_pt_result7 = np.array(mito_pt_result7[1:], dtype=float)
    mito_gc_content = np.array(mito_gc_content[1:], dtype=float)
    mito_nm_length = np.array(mito_nm_length[1:], dtype=float)
    mito_pt_length = np.array(mito_pt_length[1:], dtype=float)

    return mito_nm_result1, mito_pt_result1, mito_nm_result3, mito_pt_result3, mito_nm_result7, mito_pt_result7, mito_gc_content, mito_nm_length, mito_pt_length

def plastidData(plastid_nm_result1, plastid_pt_result1, plastid_nm_result3, plastid_pt_result3, plastid_nm_result7, plastid_pt_result7, plastid_gc_content, plastid_nm_length, plastid_pt_length):

    countNM, countPT, countGC, countLen = 0, 0, 0, 0

    plastidFileNM = dir_path + '/../Analysis/GeCo/Plastid/geco3_Plastid.csv'
    plastidFilePT = dir_path + '/../Analysis/AC/Plastid/ac2_Plastid.csv'
    plastidFileGC = dir_path + '/../Analysis/GCcontent/Plastid/gc_content_Plastid.csv'
    plastidFileLN = dir_path + '/../Analysis/LengthSeq/Plastid/length_Plastid.csv'

    if trainDatabase > 0:
        numFilesNM = int(len(open(plastidFileNM).readlines(  )) * trainDatabase)
        numFilesPT = int(len(open(plastidFilePT).readlines(  )) * trainDatabase)
        numFilesGC = int(len(open(plastidFileGC).readlines(  )) * trainDatabase)
        numFilesLen = int(len(open(plastidFileLN).readlines(  )) * trainDatabase)
    else:
        numFilesNM_ForTestMin = int(len(open(plastidFileNM).readlines(  )) * crossValidation_Min_Lim)
        numFilesNM_ForTestMax = int(len(open(plastidFileNM).readlines(  )) * crossValidation_Max_Lim)
        numFilesPT_ForTestMin = int(len(open(plastidFilePT).readlines(  )) * crossValidation_Min_Lim)
        numFilesPT_ForTestMax = int(len(open(plastidFilePT).readlines(  )) * crossValidation_Max_Lim)
        numFilesGC_ForTestMin = int(len(open(plastidFileGC).readlines(  )) * crossValidation_Min_Lim)
        numFilesGC_ForTestMax = int(len(open(plastidFileGC).readlines(  )) * crossValidation_Max_Lim)
        numFilesLN_ForTestMin = int(len(open(plastidFileLN).readlines(  )) * crossValidation_Min_Lim)
        numFilesLN_ForTestMax = int(len(open(plastidFileLN).readlines(  )) * crossValidation_Max_Lim)

    # Plastid CSV files
    with open('Analysis/GeCo/Plastid/geco3_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countNM == numFilesNM:
                    break
                else:
                    plastid_nm_result1.append("0"+row[0].split("\t")[L1])
                    plastid_nm_result3.append("0"+row[0].split("\t")[L3])
                    plastid_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1
            else:
                if countNM >= numFilesNM_ForTestMin and countNM <= numFilesNM_ForTestMax:
                    countNM = countNM + 1
                else:
                    plastid_nm_result1.append("0"+row[0].split("\t")[L1])
                    plastid_nm_result3.append("0"+row[0].split("\t")[L3])
                    plastid_nm_result7.append("0"+row[0].split("\t")[L7])
                    countNM = countNM + 1

    with open('Analysis/AC/Plastid/ac2_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countPT == numFilesPT:
                    break
                else:
                    plastid_pt_result1.append(row[0].split("\t")[L1])
                    plastid_pt_result3.append(row[0].split("\t")[L3])
                    plastid_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1
            else:
                if countPT >= numFilesPT_ForTestMin and countPT <= numFilesPT_ForTestMax:
                    countPT = countPT + 1
                else:
                    plastid_pt_result1.append(row[0].split("\t")[L1])
                    plastid_pt_result3.append(row[0].split("\t")[L3])
                    plastid_pt_result7.append(row[0].split("\t")[L7])
                    countPT = countPT + 1

    with open('Analysis/GCcontent/Plastid/gc_content_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countGC == numFilesGC:
                    break
                else:
                    plastid_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1
            else:
                if countGC >= numFilesGC_ForTestMin and countGC <= numFilesGC_ForTestMax:
                    countGC = countGC + 1
                else:
                    plastid_gc_content.append("0"+row[0].split("\t")[3])
                    countGC = countGC + 1

    with open('Analysis/LengthSeq/Plastid/length_Plastid.csv', 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if trainDatabase > 0:
                if countLen == numFilesLen:
                    break
                else:
                    plastid_nm_length.append(row[0].split("\t")[1])
                    plastid_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1
            else:
                if countLen >= numFilesLN_ForTestMin and countLen <= numFilesLN_ForTestMax:
                    countLen = countLen + 1
                else:
                    plastid_nm_length.append(row[0].split("\t")[1])
                    plastid_pt_length.append(row[0].split("\t")[2])
                    countLen = countLen + 1

    # Data conversion and header removal
    plastid_nm_result1 = np.array(plastid_nm_result1[1:], dtype=float)
    plastid_pt_result1 = np.array(plastid_pt_result1[1:], dtype=float)
    plastid_nm_result3 = np.array(plastid_nm_result3[1:], dtype=float)
    plastid_pt_result3 = np.array(plastid_pt_result3[1:], dtype=float)
    plastid_nm_result7 = np.array(plastid_nm_result7[1:], dtype=float)
    plastid_pt_result7 = np.array(plastid_pt_result7[1:], dtype=float)
    plastid_gc_content = np.array(plastid_gc_content[1:], dtype=float)
    plastid_nm_length = np.array(plastid_nm_length[1:], dtype=float)
    plastid_pt_length = np.array(plastid_pt_length[1:], dtype=float)

    return plastid_nm_result1, plastid_pt_result1, plastid_nm_result3, plastid_pt_result3, plastid_nm_result7, plastid_pt_result7, plastid_gc_content, plastid_nm_length, plastid_pt_length

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

def calcProb(p_type, likelihood_dna_1, likelihood_aa1, likelihood_dna_3, likelihood_aa3, likelihood_dna_7, likelihood_aa7, likelihood_gc, likelihood_dna_length, likelihood_aa_length):
    # TODO: Gradiente Descendente
    if preditors_required == "1111": # All predictors
        return p_type + likelihood_dna_1 + likelihood_aa1 + likelihood_dna_3 + likelihood_aa3 + likelihood_dna_7 + likelihood_aa7 + likelihood_gc + likelihood_dna_length + likelihood_aa_length
    elif preditors_required == "0001": # Only DNA Compression
        return p_type + likelihood_dna_1 + likelihood_dna_3 + likelihood_dna_7
    elif preditors_required == "0010": # Only AA Compression
        return p_type + likelihood_aa1 + likelihood_aa3 + likelihood_aa7
    elif preditors_required == "0011": # Only GC-Content
        return p_type + likelihood_gc
    elif preditors_required == "0100": # Only DNA Length
        return p_type + likelihood_dna_length
    elif preditors_required == "0101": # Only AA Length
        return p_type + likelihood_aa_length
    elif preditors_required == "0110": # DNA & AA Compression
        return p_type + likelihood_dna_1 + likelihood_aa1 + likelihood_dna_3 + likelihood_aa3 + likelihood_dna_7 + likelihood_aa7
    elif preditors_required == "0111": # DNA & AA Compression + GC-Content
        return p_type + likelihood_dna_1 + likelihood_aa1 + likelihood_dna_3 + likelihood_aa3 + likelihood_dna_7 + likelihood_aa7 + likelihood_gc
    elif preditors_required == "1000": # DNA & AA Compression + GC-Content + DNA Length
        return p_type + likelihood_dna_1 + likelihood_aa1 + likelihood_dna_3 + likelihood_aa3 + likelihood_dna_7 + likelihood_aa7 + likelihood_gc + likelihood_dna_length
    elif preditors_required == "1001": # DNA & AA Compression + GC-Content + DNA & AA Length
        return p_type + likelihood_dna_1 + likelihood_aa1 + likelihood_dna_3 + likelihood_aa3 + likelihood_dna_7 + likelihood_aa7 + likelihood_gc + likelihood_dna_length + likelihood_aa_length
    elif preditors_required == "1010": # DNA & AA Compression + DNA & AA Length
        return p_type + likelihood_dna_1 + likelihood_aa1 + likelihood_dna_3 + likelihood_aa3 + likelihood_dna_7 + likelihood_aa7 + likelihood_dna_length + likelihood_aa_length
    else:
        print("Insert a correct code for the predictors wanted!")
        print(" 1111 -> All predictors")
        print(" 0001 -> Only DNA Compression")
        print(" 0010 -> Only AA Compression")
        print(" 0011 -> Only GC-Content")
        print(" 0100 -> Only DNA Length")
        print(" 0101 -> Only AA Length")
        print(" 0110 -> DNA & AA Compression")
        print(" 0111 -> DNA & AA Compression + GC-Content")
        print(" 1000 -> DNA & AA Compression + GC-Content + DNA Length")
        print(" 1001 -> DNA & AA Compression + GC-Content + DNA & AA Length")
        print(" 1010 -> DNA & AA Compression + DNA & AA Length")
        return 0;

def domainAnalysis(probabilities):
    domains=["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]
    return domains[probabilities.index(max(probabilities))]

virus_nm_result1, virus_pt_result1, virus_nm_result3, virus_pt_result3, virus_nm_result7, virus_pt_result7, virus_gc_content, virus_nm_length, virus_pt_length = virusData(virus_nm_result1, virus_pt_result1, virus_nm_result3, virus_pt_result3, virus_nm_result7, virus_pt_result7, virus_gc_content, virus_nm_length, virus_pt_length)
bacteria_nm_result1, bacteria_pt_result1, bacteria_nm_result3, bacteria_pt_result3, bacteria_nm_result7, bacteria_pt_result7, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length = bacteriaData(bacteria_nm_result1, bacteria_pt_result1, bacteria_nm_result3, bacteria_pt_result3, bacteria_nm_result7, bacteria_pt_result7, bacteria_gc_content, bacteria_nm_length, bacteria_pt_length)
archaea_nm_result1, archaea_pt_result1, archaea_nm_result3, archaea_pt_result3, archaea_nm_result7, archaea_pt_result7, archaea_gc_content, archaea_nm_length, archaea_pt_length = archaeaData(archaea_nm_result1, archaea_pt_result1, archaea_nm_result3, archaea_pt_result3, archaea_nm_result7, archaea_pt_result7, archaea_gc_content, archaea_nm_length, archaea_pt_length)
fungi_nm_result1, fungi_pt_result1, fungi_nm_result3, fungi_pt_result3, fungi_nm_result7, fungi_pt_result7, fungi_gc_content, fungi_nm_length, fungi_pt_length = fungiData(fungi_nm_result1, fungi_pt_result1, fungi_nm_result3, fungi_pt_result3, fungi_nm_result7, fungi_pt_result7, fungi_gc_content, fungi_nm_length, fungi_pt_length)
plant_nm_result1, plant_pt_result1, plant_nm_result3, plant_pt_result3, plant_nm_result7, plant_pt_result7, plant_gc_content, plant_nm_length, plant_pt_length = plantData(plant_nm_result1, plant_pt_result1, plant_nm_result3, plant_pt_result3, plant_nm_result7, plant_pt_result7, plant_gc_content, plant_nm_length, plant_pt_length)
protozoa_nm_result1, protozoa_pt_result1, protozoa_nm_result3, protozoa_pt_result3, protozoa_nm_result7, protozoa_pt_result7, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length = protozoaData(protozoa_nm_result1, protozoa_pt_result1, protozoa_nm_result3, protozoa_pt_result3, protozoa_nm_result7, protozoa_pt_result7, protozoa_gc_content, protozoa_nm_length, protozoa_pt_length)
mito_nm_result1, mito_pt_result1, mito_nm_result3, mito_pt_result3, mito_nm_result7, mito_pt_result7, mito_gc_content, mito_nm_length, mito_pt_length = mitoData(mito_nm_result1, mito_pt_result1, mito_nm_result3, mito_pt_result3, mito_nm_result7, mito_pt_result7, mito_gc_content, mito_nm_length, mito_pt_length)
plastid_nm_result1, plastid_pt_result1, plastid_nm_result3, plastid_pt_result3, plastid_nm_result7, plastid_pt_result7, plastid_gc_content, plastid_nm_length, plastid_pt_length = plastidData(plastid_nm_result1, plastid_pt_result1, plastid_nm_result3, plastid_pt_result3, plastid_nm_result7, plastid_pt_result7, plastid_gc_content, plastid_nm_length, plastid_pt_length)

virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length, maxvaluefound_length_DNA = normalizeLengths(virus_nm_length, bacteria_nm_length, archaea_nm_length, fungi_nm_length, plant_nm_length, protozoa_nm_length, mito_nm_length, plastid_nm_length)
virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length, maxvaluefound_length_AA = normalizeLengths(virus_pt_length, bacteria_pt_length, archaea_pt_length, fungi_pt_length, plant_pt_length, protozoa_pt_length, mito_pt_length, plastid_pt_length)

# Fit a normal distribution to the data
mu_virus_nm1, std_virus_nm1 = norm.fit(virus_nm_result1)
mu_virus_pt1, std_virus_pt1 = norm.fit(virus_pt_result1)
mu_virus_nm3, std_virus_nm3 = norm.fit(virus_nm_result3)
mu_virus_pt3, std_virus_pt3 = norm.fit(virus_pt_result3)
mu_virus_nm7, std_virus_nm7 = norm.fit(virus_nm_result7)
mu_virus_pt7, std_virus_pt7 = norm.fit(virus_pt_result7)
mu_virus_gc, std_virus_gc = norm.fit(virus_gc_content)
mu_virus_nmL, std_virus_nmL = norm.fit(virus_nm_length)
mu_virus_ptL, std_virus_ptL = norm.fit(virus_pt_length)

mu_bacteria_nm1, std_bacteria_nm1 = norm.fit(bacteria_nm_result1)
mu_bacteria_pt1, std_bacteria_pt1 = norm.fit(bacteria_pt_result1)
mu_bacteria_nm3, std_bacteria_nm3 = norm.fit(bacteria_nm_result3)
mu_bacteria_pt3, std_bacteria_pt3 = norm.fit(bacteria_pt_result3)
mu_bacteria_nm7, std_bacteria_nm7 = norm.fit(bacteria_nm_result7)
mu_bacteria_pt7, std_bacteria_pt7 = norm.fit(bacteria_pt_result7)
mu_bacteria_gc, std_bacteria_gc = norm.fit(bacteria_gc_content)
mu_bacteria_nmL, std_bacteria_nmL = norm.fit(bacteria_nm_length)
mu_bacteria_ptL, std_bacteria_ptL = norm.fit(bacteria_pt_length)

mu_archaea_nm1, std_archaea_nm1 = norm.fit(archaea_nm_result1)
mu_archaea_pt1, std_archaea_pt1 = norm.fit(archaea_pt_result1)
mu_archaea_nm3, std_archaea_nm3 = norm.fit(archaea_nm_result3)
mu_archaea_pt3, std_archaea_pt3 = norm.fit(archaea_pt_result3)
mu_archaea_nm7, std_archaea_nm7 = norm.fit(archaea_nm_result7)
mu_archaea_pt7, std_archaea_pt7 = norm.fit(archaea_pt_result7)
mu_archaea_gc, std_archaea_gc = norm.fit(archaea_gc_content)
mu_archaea_nmL, std_archaea_nmL = norm.fit(archaea_nm_length)
mu_archaea_ptL, std_archaea_ptL = norm.fit(archaea_pt_length)

mu_fungi_nm1, std_fungi_nm1 = norm.fit(fungi_nm_result1)
mu_fungi_pt1, std_fungi_pt1 = norm.fit(fungi_pt_result1)
mu_fungi_nm3, std_fungi_nm3 = norm.fit(fungi_nm_result3)
mu_fungi_pt3, std_fungi_pt3 = norm.fit(fungi_pt_result3)
mu_fungi_nm7, std_fungi_nm7 = norm.fit(fungi_nm_result7)
mu_fungi_pt7, std_fungi_pt7 = norm.fit(fungi_pt_result7)
mu_fungi_gc, std_fungi_gc = norm.fit(fungi_gc_content)
mu_fungi_nmL, std_fungi_nmL = norm.fit(fungi_nm_length)
mu_fungi_ptL, std_fungi_ptL = norm.fit(fungi_pt_length)

mu_plant_nm1, std_plant_nm1 = norm.fit(plant_nm_result1)
mu_plant_pt1, std_plant_pt1 = norm.fit(plant_pt_result1)
mu_plant_nm3, std_plant_nm3 = norm.fit(plant_nm_result3)
mu_plant_pt3, std_plant_pt3 = norm.fit(plant_pt_result3)
mu_plant_nm7, std_plant_nm7 = norm.fit(plant_nm_result7)
mu_plant_pt7, std_plant_pt7 = norm.fit(plant_pt_result7)
mu_plant_gc, std_plant_gc = norm.fit(plant_gc_content)
mu_plant_nmL, std_plant_nmL = norm.fit(plant_nm_length)
mu_plant_ptL, std_plant_ptL = norm.fit(plant_pt_length)

mu_protozoa_nm1, std_protozoa_nm1 = norm.fit(protozoa_nm_result1)
mu_protozoa_pt1, std_protozoa_pt1 = norm.fit(protozoa_pt_result1)
mu_protozoa_nm3, std_protozoa_nm3 = norm.fit(protozoa_nm_result3)
mu_protozoa_pt3, std_protozoa_pt3 = norm.fit(protozoa_pt_result3)
mu_protozoa_nm7, std_protozoa_nm7 = norm.fit(protozoa_nm_result7)
mu_protozoa_pt7, std_protozoa_pt7 = norm.fit(protozoa_pt_result7)
mu_protozoa_gc, std_protozoa_gc = norm.fit(protozoa_gc_content)
mu_protozoa_nmL, std_protozoa_nmL = norm.fit(protozoa_nm_length)
mu_protozoa_ptL, std_protozoa_ptL = norm.fit(protozoa_pt_length)

mu_mito_nm1, std_mito_nm1 = norm.fit(mito_nm_result1)
mu_mito_pt1, std_mito_pt1 = norm.fit(mito_pt_result1)
mu_mito_nm3, std_mito_nm3 = norm.fit(mito_nm_result3)
mu_mito_pt3, std_mito_pt3 = norm.fit(mito_pt_result3)
mu_mito_nm7, std_mito_nm7 = norm.fit(mito_nm_result7)
mu_mito_pt7, std_mito_pt7 = norm.fit(mito_pt_result7)
mu_mito_gc, std_mito_gc = norm.fit(mito_gc_content)
mu_mito_nmL, std_mito_nmL = norm.fit(mito_nm_length)
mu_mito_ptL, std_mito_ptL = norm.fit(mito_pt_length)

mu_plastid_nm1, std_plastid_nm1 = norm.fit(plastid_nm_result1)
mu_plastid_pt1, std_plastid_pt1 = norm.fit(plastid_pt_result1)
mu_plastid_nm3, std_plastid_nm3 = norm.fit(plastid_nm_result3)
mu_plastid_pt3, std_plastid_pt3 = norm.fit(plastid_pt_result3)
mu_plastid_nm7, std_plastid_nm7 = norm.fit(plastid_nm_result7)
mu_plastid_pt7, std_plastid_pt7 = norm.fit(plastid_pt_result7)
mu_plastid_gc, std_plastid_gc = norm.fit(plastid_gc_content)
mu_plastid_nmL, std_plastid_nmL = norm.fit(plastid_nm_length)
mu_plastid_ptL, std_plastid_ptL = norm.fit(plastid_pt_length)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                    Probabilities (Log base 2)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

p_type = math.log(1/N,2)

likelihood_virus_dna1 = math.log(NormalDist(mu_virus_nm1,std_virus_nm1).pdf(sample_DNA_1),2)
likelihood_virus_aa1 = math.log(NormalDist(mu_virus_pt1,std_virus_pt1).pdf(sample_AA_1),2)
likelihood_virus_dna3 = math.log(NormalDist(mu_virus_nm3,std_virus_nm3).pdf(sample_DNA_3),2)
likelihood_virus_aa3 = math.log(NormalDist(mu_virus_pt3,std_virus_pt3).pdf(sample_AA_3),2)
likelihood_virus_dna7 = math.log(NormalDist(mu_virus_nm7,std_virus_nm7).pdf(sample_DNA_7),2)
likelihood_virus_aa7 = math.log(NormalDist(mu_virus_pt7,std_virus_pt7).pdf(sample_AA_7),2)
likelihood_virus_gc = math.log(NormalDist(mu_virus_gc,std_virus_gc).pdf(gc_percent),2)
likelihood_virus_dnaL = math.log(NormalDist(mu_virus_nmL,std_virus_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_virus_aaL = math.log(NormalDist(mu_virus_ptL,std_virus_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_bacteria_dna1 = math.log(NormalDist(mu_bacteria_nm1,std_bacteria_nm1).pdf(sample_DNA_1),2)
likelihood_bacteria_aa1 = math.log(NormalDist(mu_bacteria_pt1,std_bacteria_pt1).pdf(sample_AA_1),2)
likelihood_bacteria_dna3 = math.log(NormalDist(mu_bacteria_nm3,std_bacteria_nm3).pdf(sample_DNA_3),2)
likelihood_bacteria_aa3 = math.log(NormalDist(mu_bacteria_pt3,std_bacteria_pt3).pdf(sample_AA_3),2)
likelihood_bacteria_dna7 = math.log(NormalDist(mu_bacteria_nm7,std_bacteria_nm7).pdf(sample_DNA_7),2)
likelihood_bacteria_aa7 = math.log(NormalDist(mu_bacteria_pt7,std_bacteria_pt7).pdf(sample_AA_7),2)
likelihood_bacteria_gc = math.log(NormalDist(mu_bacteria_gc,std_bacteria_gc).pdf(gc_percent),2)
likelihood_bacteria_dnaL = math.log(NormalDist(mu_bacteria_nmL,std_bacteria_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_bacteria_aaL = math.log(NormalDist(mu_bacteria_ptL,std_bacteria_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_archaea_dna1 = math.log(NormalDist(mu_archaea_nm1,std_archaea_nm1).pdf(sample_DNA_1),2)
likelihood_archaea_aa1 = math.log(NormalDist(mu_archaea_pt1,std_archaea_pt1).pdf(sample_AA_1),2)
likelihood_archaea_dna3 = math.log(NormalDist(mu_archaea_nm3,std_archaea_nm3).pdf(sample_DNA_3),2)
likelihood_archaea_aa3 = math.log(NormalDist(mu_archaea_pt3,std_archaea_pt3).pdf(sample_AA_3),2)
likelihood_archaea_dna7 = math.log(NormalDist(mu_archaea_nm7,std_archaea_nm7).pdf(sample_DNA_7),2)
likelihood_archaea_aa7 = math.log(NormalDist(mu_archaea_pt7,std_archaea_pt7).pdf(sample_AA_7),2)
likelihood_archaea_gc = math.log(NormalDist(mu_archaea_gc,std_archaea_gc).pdf(gc_percent),2)
likelihood_archaea_dnaL = math.log(NormalDist(mu_archaea_nmL,std_archaea_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_archaea_aaL = math.log(NormalDist(mu_archaea_ptL,std_archaea_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_fungi_dna1 = math.log(NormalDist(mu_fungi_nm1,std_fungi_nm1).pdf(sample_DNA_1),2)
likelihood_fungi_aa1 = math.log(NormalDist(mu_fungi_pt1,std_fungi_pt1).pdf(sample_AA_1),2)
likelihood_fungi_dna3 = math.log(NormalDist(mu_fungi_nm3,std_fungi_nm3).pdf(sample_DNA_3),2)
likelihood_fungi_aa3 = math.log(NormalDist(mu_fungi_pt3,std_fungi_pt3).pdf(sample_AA_3),2)
likelihood_fungi_dna7 = math.log(NormalDist(mu_fungi_nm7,std_fungi_nm7).pdf(sample_DNA_7),2)
likelihood_fungi_aa7 = math.log(NormalDist(mu_fungi_pt7,std_fungi_pt7).pdf(sample_AA_7),2)
likelihood_fungi_gc = math.log(NormalDist(mu_fungi_gc,std_fungi_gc).pdf(gc_percent),2)
likelihood_fungi_dnaL = math.log(NormalDist(mu_fungi_nmL,std_fungi_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_fungi_aaL = math.log(NormalDist(mu_fungi_ptL,std_fungi_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_plant_dna1 = math.log(NormalDist(mu_plant_nm1,std_plant_nm1).pdf(sample_DNA_1),2)
likelihood_plant_aa1 = math.log(NormalDist(mu_plant_pt1,std_plant_pt1).pdf(sample_AA_1),2)
likelihood_plant_dna3 = math.log(NormalDist(mu_plant_nm3,std_plant_nm3).pdf(sample_DNA_3),2)
likelihood_plant_aa3 = math.log(NormalDist(mu_plant_pt3,std_plant_pt3).pdf(sample_AA_3),2)
likelihood_plant_dna7 = math.log(NormalDist(mu_plant_nm7,std_plant_nm7).pdf(sample_DNA_7),2)
likelihood_plant_aa7 = math.log(NormalDist(mu_plant_pt7,std_plant_pt7).pdf(sample_AA_7),2)
likelihood_plant_gc = math.log(NormalDist(mu_plant_gc,std_plant_gc).pdf(gc_percent),2)
likelihood_plant_dnaL = math.log(NormalDist(mu_plant_nmL,std_plant_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_plant_aaL = math.log(NormalDist(mu_plant_ptL,std_plant_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_protozoa_dna1 = math.log(NormalDist(mu_protozoa_nm1,std_protozoa_nm1).pdf(sample_DNA_1),2)
likelihood_protozoa_aa1 = math.log(NormalDist(mu_protozoa_pt1,std_protozoa_pt1).pdf(sample_AA_1),2)
likelihood_protozoa_dna3 = math.log(NormalDist(mu_protozoa_nm3,std_protozoa_nm3).pdf(sample_DNA_3),2)
likelihood_protozoa_aa3 = math.log(NormalDist(mu_protozoa_pt3,std_protozoa_pt3).pdf(sample_AA_3),2)
likelihood_protozoa_dna7 = math.log(NormalDist(mu_protozoa_nm7,std_protozoa_nm7).pdf(sample_DNA_7),2)
likelihood_protozoa_aa7 = math.log(NormalDist(mu_protozoa_pt7,std_protozoa_pt7).pdf(sample_AA_7),2)
likelihood_protozoa_gc = math.log(NormalDist(mu_protozoa_gc,std_protozoa_gc).pdf(gc_percent),2)
likelihood_protozoa_dnaL = math.log(NormalDist(mu_protozoa_nmL,std_protozoa_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_protozoa_aaL = math.log(NormalDist(mu_protozoa_ptL,std_protozoa_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_mito_dna1 = math.log(NormalDist(mu_mito_nm1,std_mito_nm1).pdf(sample_DNA_1),2)
likelihood_mito_aa1 = math.log(NormalDist(mu_mito_pt1,std_mito_pt1).pdf(sample_AA_1),2)
likelihood_mito_dna3 = math.log(NormalDist(mu_mito_nm3,std_mito_nm3).pdf(sample_DNA_3),2)
likelihood_mito_aa3 = math.log(NormalDist(mu_mito_pt3,std_mito_pt3).pdf(sample_AA_3),2)
likelihood_mito_dna7 = math.log(NormalDist(mu_mito_nm7,std_mito_nm7).pdf(sample_DNA_7),2)
likelihood_mito_aa7 = math.log(NormalDist(mu_mito_pt7,std_mito_pt7).pdf(sample_AA_7),2)
likelihood_mito_gc = math.log(NormalDist(mu_mito_gc,std_mito_gc).pdf(gc_percent),2)
likelihood_mito_dnaL = math.log(NormalDist(mu_mito_nmL,std_mito_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_mito_aaL = math.log(NormalDist(mu_mito_ptL,std_mito_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

likelihood_plastid_dna1 = math.log(NormalDist(mu_plastid_nm1,std_plastid_nm1).pdf(sample_DNA_1),2)
likelihood_plastid_aa1 = math.log(NormalDist(mu_plastid_pt1,std_plastid_pt1).pdf(sample_AA_1),2)
likelihood_plastid_dna3 = math.log(NormalDist(mu_plastid_nm3,std_plastid_nm3).pdf(sample_DNA_3),2)
likelihood_plastid_aa3 = math.log(NormalDist(mu_plastid_pt3,std_plastid_pt3).pdf(sample_AA_3),2)
likelihood_plastid_dna7 = math.log(NormalDist(mu_plastid_nm7,std_plastid_nm7).pdf(sample_DNA_7),2)
likelihood_plastid_aa7 = math.log(NormalDist(mu_plastid_pt7,std_plastid_pt7).pdf(sample_AA_7),2)
likelihood_plastid_gc = math.log(NormalDist(mu_plastid_gc,std_plastid_gc).pdf(gc_percent),2)
likelihood_plastid_dnaL = math.log(NormalDist(mu_plastid_nmL,std_plastid_nmL).pdf(float("{:.10f}".format(sample_DNA_LEN/maxvaluefound_length_DNA)))+1,2)
likelihood_plastid_aaL = math.log(NormalDist(mu_plastid_ptL,std_plastid_ptL).pdf(float("{:.10f}".format(sample_AA_LEN/maxvaluefound_length_AA)))+1,2)

# Final Probabilities
p_virus = calcProb(p_type, likelihood_virus_dna1, likelihood_virus_aa1, likelihood_virus_dna3, likelihood_virus_aa3, likelihood_virus_dna7, likelihood_virus_aa7, likelihood_virus_gc, likelihood_virus_dnaL, likelihood_virus_aaL)
p_bacteria = calcProb(p_type, likelihood_bacteria_dna1, likelihood_bacteria_aa1, likelihood_bacteria_dna3, likelihood_bacteria_aa3, likelihood_bacteria_dna7, likelihood_bacteria_aa7, likelihood_bacteria_gc, likelihood_bacteria_dnaL, likelihood_bacteria_aaL)
p_archaea = calcProb(p_type, likelihood_archaea_dna1, likelihood_archaea_aa1, likelihood_archaea_dna3, likelihood_archaea_aa3, likelihood_archaea_dna7, likelihood_archaea_aa7, likelihood_archaea_gc, likelihood_archaea_dnaL, likelihood_archaea_aaL)
p_fungi = calcProb(p_type, likelihood_fungi_dna1, likelihood_fungi_aa1, likelihood_fungi_dna3, likelihood_fungi_aa3, likelihood_fungi_dna7, likelihood_fungi_aa7, likelihood_fungi_gc, likelihood_fungi_dnaL, likelihood_fungi_aaL)
p_plant = calcProb(p_type, likelihood_plant_dna1, likelihood_plant_aa1, likelihood_plant_dna3, likelihood_plant_aa3, likelihood_plant_dna7, likelihood_plant_aa7, likelihood_plant_gc, likelihood_plant_dnaL, likelihood_plant_aaL)
p_protozoa = calcProb(p_type, likelihood_protozoa_dna1, likelihood_protozoa_aa1, likelihood_protozoa_dna3, likelihood_protozoa_aa3, likelihood_protozoa_dna7, likelihood_protozoa_aa7, likelihood_protozoa_gc, likelihood_protozoa_dnaL, likelihood_protozoa_aaL)
p_mito = calcProb(p_type, likelihood_mito_dna1, likelihood_mito_aa1, likelihood_mito_dna3, likelihood_mito_aa3, likelihood_mito_dna7, likelihood_mito_aa7, likelihood_mito_gc, likelihood_mito_dnaL, likelihood_mito_aaL)
p_plastid = calcProb(p_type, likelihood_plastid_dna1, likelihood_plastid_aa1, likelihood_plastid_dna3, likelihood_plastid_aa3, likelihood_plastid_dna7, likelihood_plastid_aa7, likelihood_plastid_gc, likelihood_plastid_dnaL, likelihood_plastid_aaL)

probabilities=[p_virus, p_bacteria, p_archaea, p_fungi, p_plant, p_protozoa, p_mito, p_plastid]

#print(probabilities)
print(domainAnalysis(probabilities))

#generateHist_GC()

#generateHist(virus_nm_result3,mu_virus_nm3,std_virus_nm3)