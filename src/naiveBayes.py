import sys, csv, math
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
from statistics import NormalDist

# ( https://www.youtube.com/watch?v=rzFX5NWojp0 ) 

N = 8                               # Number of types (Domains)
sample_DNA = float(sys.argv[1])     # Normalize Compression Rate of the DNA for the Input Sequence
sample_AA = float(sys.argv[2])      # Normalize Compression Rate of the AA for the Input Sequence
gc_percent = float(sys.argv[3])     # GC-Content Percentage ([0..1]) of the Input Sequence

## Virus Values
virus_nm_result7, virus_pt_result7, virus_gc_content = [], [], []

## Bacterias Values
bacteria_nm_result7, bacteria_pt_result7, bacteria_gc_content = [], [], []

## Archaeas Values
archaea_nm_result7, archaea_pt_result7, archaea_gc_content = [], [], []

## Fungis Values
fungi_nm_result7, fungi_pt_result7, fungi_gc_content = [], [], []

## Plants Values
plant_nm_result7, plant_pt_result7, plant_gc_content = [], [], []

## Protozoas Values
protozoa_nm_result7, protozoa_pt_result7, protozoa_gc_content = [], [], []

## Mitochondrials Values
mito_nm_result7, mito_pt_result7, mito_gc_content = [], [], []

## Platids Values
plastid_nm_result7, plastid_pt_result7, plastid_gc_content = [], [], []

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

def domainAnalysis(probabilities):
    domains=["Virus", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]
    return domains[probabilities.index(max(probabilities))]

# Viral CSV files
with open('GeCo3_Output/Virus/NM/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Virus/PT/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Virus/gc_content_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_gc_content.append("0"+row[0].split("\t")[3])

'''
# Bacterias CSV files
with open('GeCo3_Output/Bacteria/NM/geco3_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacteria_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Bacteria/PT/geco3_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacteria_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Bacteria/gc_content_Bacteria.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        bacteria_gc_content.append("0"+row[0].split("\t")[3])
'''

# Archaeas CSV files
with open('GeCo3_Output/Archaea/NM/geco3_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaea_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Archaea/PT/geco3_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaea_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Archaea/gc_content_Archaea.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        archaea_gc_content.append("0"+row[0].split("\t")[3])

# Fungis CSV files
with open('GeCo3_Output/Fungi/NM/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Fungi/PT/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Fungi/gc_content_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_gc_content.append("0"+row[0].split("\t")[3])

# Plants CSV files
with open('GeCo3_Output/Plant/NM/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Plant/PT/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Plant/gc_content_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_gc_content.append("0"+row[0].split("\t")[3])

# Protozoa CSV files
with open('GeCo3_Output/Protozoa/NM/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Protozoa/PT/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Protozoa/gc_content_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_gc_content.append("0"+row[0].split("\t")[3])

# Mitochondrial CSV files
with open('GeCo3_Output/Mitochondrial/NM/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Mitochondrial/PT/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Mitochondrial/gc_content_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_gc_content.append("0"+row[0].split("\t")[3])

# Plastid CSV files
with open('GeCo3_Output/Plastid/NM/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Plastid/PT/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_pt_result7.append("0"+row[0].split("\t")[12])

with open('Analysis/GCcontent/Plastid/gc_content_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_gc_content.append("0"+row[0].split("\t")[3])

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 
# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

# Data conversion and header removal
virus_nm_result7 = np.array(virus_nm_result7[1:], dtype=float)
virus_pt_result7 = np.array(virus_pt_result7[1:], dtype=float)
virus_gc_content = np.array(virus_gc_content[1:], dtype=float)

#bacteria_nm_result7 = np.array(bacteria_nm_result7[1:], dtype=float)
#bacteria_pt_result7 = np.array(bacteria_pt_result7[1:], dtype=float)
#bacteria_gc_content = np.array(bacteria_gc_content[1:], dtype=float)

archaea_nm_result7 = np.array(archaea_nm_result7[1:], dtype=float)
archaea_pt_result7 = np.array(archaea_pt_result7[1:], dtype=float)
archaea_gc_content = np.array(archaea_gc_content[1:], dtype=float)

fungi_nm_result7 = np.array(fungi_nm_result7[1:], dtype=float)
fungi_pt_result7 = np.array(fungi_pt_result7[1:], dtype=float)
fungi_gc_content = np.array(fungi_gc_content[1:], dtype=float)

plant_nm_result7 = np.array(plant_nm_result7[1:], dtype=float)
plant_pt_result7 = np.array(plant_pt_result7[1:], dtype=float)
plant_gc_content = np.array(plant_gc_content[1:], dtype=float)

protozoa_nm_result7 = np.array(protozoa_nm_result7[1:], dtype=float)
protozoa_pt_result7 = np.array(protozoa_pt_result7[1:], dtype=float)
protozoa_gc_content = np.array(protozoa_gc_content[1:], dtype=float)

mito_nm_result7 = np.array(mito_nm_result7[1:], dtype=float)
mito_pt_result7 = np.array(mito_pt_result7[1:], dtype=float)
mito_gc_content = np.array(mito_gc_content[1:], dtype=float)

plastid_nm_result7 = np.array(plastid_nm_result7[1:], dtype=float)
plastid_pt_result7 = np.array(plastid_pt_result7[1:], dtype=float)
plastid_gc_content = np.array(plastid_gc_content[1:], dtype=float)

# Fit a normal distribution to the data
mu_virus_nm, std_virus_nm = norm.fit(virus_nm_result7)
mu_virus_pt, std_virus_pt = norm.fit(virus_pt_result7)
mu_virus_gc, std_virus_gc = norm.fit(virus_gc_content)

#mu_bacteria_nm, std_bacteria_nm = norm.fit(bacteria_nm_result7)
#mu_bacteria_pt, std_bacteria_pt = norm.fit(bacteria_pt_result7)
#mu_bacteria_gc, std_bacteria_gc = norm.fit(bacteria_gc_content)

mu_archaea_nm, std_archaea_nm = norm.fit(archaea_nm_result7)
mu_archaea_pt, std_archaea_pt = norm.fit(archaea_pt_result7)
mu_archaea_gc, std_archaea_gc = norm.fit(archaea_gc_content)

mu_fungi_nm, std_fungi_nm = norm.fit(fungi_nm_result7)
mu_fungi_pt, std_fungi_pt = norm.fit(fungi_pt_result7)
mu_fungi_gc, std_fungi_gc = norm.fit(fungi_gc_content)

mu_plant_nm, std_plant_nm = norm.fit(plant_nm_result7)
mu_plant_pt, std_plant_pt = norm.fit(plant_pt_result7)
mu_plant_gc, std_plant_gc = norm.fit(plant_gc_content)

mu_protozoa_nm, std_protozoa_nm = norm.fit(protozoa_nm_result7)
mu_protozoa_pt, std_protozoa_pt = norm.fit(protozoa_pt_result7)
mu_protozoa_gc, std_protozoa_gc = norm.fit(protozoa_gc_content)

mu_mito_nm, std_mito_nm = norm.fit(mito_nm_result7)
mu_mito_pt, std_mito_pt = norm.fit(mito_pt_result7)
mu_mito_gc, std_mito_gc = norm.fit(mito_gc_content)

mu_plastid_nm, std_plastid_nm = norm.fit(plastid_nm_result7)
mu_plastid_pt, std_plastid_pt = norm.fit(plastid_pt_result7)
mu_plastid_gc, std_plastid_gc = norm.fit(plastid_gc_content)

'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
                    Probabilities (Log base 2)
'''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''

p_type = math.log(1/N,2)

likelihood_virus_dna = math.log(NormalDist(mu_virus_nm,std_virus_nm).pdf(sample_DNA),2)
likelihood_virus_aa = math.log(NormalDist(mu_virus_pt,std_virus_pt).pdf(sample_AA),2)
likelihood_virus_gc = math.log(NormalDist(mu_virus_gc,std_virus_gc).pdf(gc_percent),2)

#likelihood_bacteria_dna = math.log(NormalDist(mu_bacteria_nm,std_bacteria_nm).pdf(sample_DNA),2)
#likelihood_bacteria_aa = math.log(NormalDist(mu_bacteria_pt,std_bacteria_pt).pdf(sample_AA),2)
#likelihood_bacteria_gc = math.log(NormalDist(mu_bacteria_gc,std_bacteria_gc).pdf(gc_percent),2)

likelihood_archaea_dna = math.log(NormalDist(mu_archaea_nm,std_archaea_nm).pdf(sample_DNA),2)
likelihood_archaea_aa = math.log(NormalDist(mu_archaea_pt,std_archaea_pt).pdf(sample_AA),2)
likelihood_archaea_gc = math.log(NormalDist(mu_archaea_gc,std_archaea_gc).pdf(gc_percent),2)

likelihood_fungi_dna = math.log(NormalDist(mu_fungi_nm,std_fungi_nm).pdf(sample_DNA),2)
likelihood_fungi_aa = math.log(NormalDist(mu_fungi_pt,std_fungi_pt).pdf(sample_AA),2)
likelihood_fungi_gc = math.log(NormalDist(mu_fungi_gc,std_fungi_gc).pdf(gc_percent),2)

likelihood_plant_dna = math.log(NormalDist(mu_plant_nm,std_plant_nm).pdf(sample_DNA),2)
likelihood_plant_aa = math.log(NormalDist(mu_plant_pt,std_plant_pt).pdf(sample_AA),2)
likelihood_plant_gc = math.log(NormalDist(mu_plant_gc,std_plant_gc).pdf(gc_percent),2)

likelihood_protozoa_dna = math.log(NormalDist(mu_protozoa_nm,std_protozoa_nm).pdf(sample_DNA),2)
likelihood_protozoa_aa = math.log(NormalDist(mu_protozoa_pt,std_protozoa_pt).pdf(sample_AA),2)
likelihood_protozoa_gc = math.log(NormalDist(mu_protozoa_gc,std_protozoa_gc).pdf(gc_percent),2)

likelihood_mito_dna = math.log(NormalDist(mu_mito_nm,std_mito_nm).pdf(sample_DNA),2)
likelihood_mito_aa = math.log(NormalDist(mu_mito_pt,std_mito_pt).pdf(sample_AA),2)
likelihood_mito_gc = math.log(NormalDist(mu_mito_gc,std_mito_gc).pdf(gc_percent),2)

likelihood_plastid_dna = math.log(NormalDist(mu_plastid_nm,std_plastid_nm).pdf(sample_DNA),2)
likelihood_plastid_aa = math.log(NormalDist(mu_plastid_pt,std_plastid_pt).pdf(sample_AA),2)
likelihood_plastid_gc = math.log(NormalDist(mu_plastid_gc,std_plastid_gc).pdf(gc_percent),2)

# Final Probabilities
p_virus = p_type + likelihood_virus_dna + likelihood_virus_aa + likelihood_virus_gc
#p_bacteria = p_type + likelihood_bacteria_dna + likelihood_bacteria_aa + likelihood_bacteria_gc
p_archaea = p_type + likelihood_archaea_dna + likelihood_archaea_aa + likelihood_archaea_gc
p_fungi = p_type + likelihood_fungi_dna + likelihood_fungi_aa + likelihood_fungi_gc
p_plant = p_type + likelihood_plant_dna + likelihood_plant_aa + likelihood_plant_gc
p_protozoa = p_type + likelihood_protozoa_dna + likelihood_protozoa_aa + likelihood_protozoa_gc
p_mito = p_type + likelihood_mito_dna + likelihood_mito_aa + likelihood_mito_gc
p_plastid = p_type + likelihood_plastid_dna + likelihood_plastid_aa + likelihood_plastid_gc

probabilities=[p_virus, p_archaea, p_fungi, p_plant, p_protozoa, p_mito, p_plastid]

print(probabilities)
print(domainAnalysis(probabilities))