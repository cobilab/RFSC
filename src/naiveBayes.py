import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
from statistics import NormalDist

# ( https://www.youtube.com/watch?v=rzFX5NWojp0 ) 

# Number of types (Domains)
N = 8

## Virus Values
virus_nm_result7, virus_pt_result7 = [], []

## Bacterias Values
bacteria_nm_result7, bacteria_pt_result7 = [], []

## Archaeas Values
archaea_nm_result7, archaea_pt_result7 = [], []

## Fungis Values
fungi_nm_result7, fungi_pt_result7 = [], []

## Plants Values
plant_nm_result7, plant_pt_result7 = [], []

## Protozoas Values
protozoa_nm_result7, protozoa_pt_result7 = [], []

## Mitochondrials Values
mito_nm_result7, mito_pt_result7 = [], []

## Platids Values
plastid_nm_result7, plastid_pt_result7 = [], []


# Viral CSV files
with open('GeCo3_Output/Virus/NM/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Virus/PT/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_pt_result7.append("0"+row[0].split("\t")[12])

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

# Fungis CSV files
with open('GeCo3_Output/Fungi/NM/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Fungi/PT/geco3_Fungi.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        fungi_pt_result7.append("0"+row[0].split("\t")[12])

# Plants CSV files
with open('GeCo3_Output/Plant/NM/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Plant/PT/geco3_Plant.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plant_pt_result7.append("0"+row[0].split("\t")[12])

# Protozoa CSV files
with open('GeCo3_Output/Protozoa/NM/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Protozoa/PT/geco3_Protozoa.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        protozoa_pt_result7.append("0"+row[0].split("\t")[12])

# Mitochondrial CSV files
with open('GeCo3_Output/Mitochondrial/NM/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Mitochondrial/PT/geco3_Mitochondrial.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        mito_pt_result7.append("0"+row[0].split("\t")[12])

# Plastid CSV files
with open('GeCo3_Output/Plastid/NM/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_nm_result7.append("0"+row[0].split("\t")[12])

with open('GeCo3_Output/Plastid/PT/geco3_Plastid.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        plastid_pt_result7.append("0"+row[0].split("\t")[12])

# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 
# /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ /\ 
# \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ \/ 

# Generate some data for this demonstration.
#data = norm.rvs(10.0, 2.5, size=500)
#data = [8, 12, 7, 6, 2, 6, 7, 4, 3, 2, 5]
virus_nm_result7 = np.array(virus_nm_result7[1:], dtype=float)
virus_pt_result7 = np.array(virus_pt_result7[1:], dtype=float)
bacteria_nm_result7 = np.array(bacteria_nm_result7[1:], dtype=float)
bacteria_pt_result7 = np.array(bacteria_pt_result7[1:], dtype=float)
fungi_nm_result7 = np.array(fungi_nm_result7[1:], dtype=float)
fungi_pt_result7 = np.array(fungi_pt_result7[1:], dtype=float)
plant_nm_result7 = np.array(plant_nm_result7[1:], dtype=float)
plant_pt_result7 = np.array(plant_pt_result7[1:], dtype=float)
protozoa_nm_result7 = np.array(protozoa_nm_result7[1:], dtype=float)
protozoa_pt_result7 = np.array(protozoa_pt_result7[1:], dtype=float)
mito_nm_result7 = np.array(mito_nm_result7[1:], dtype=float)
mito_pt_result7 = np.array(mito_pt_result7[1:], dtype=float)
plastid_nm_result7 = np.array(plastid_nm_result7[1:], dtype=float)
plastid_pt_result7 = np.array(plastid_pt_result7[1:], dtype=float)

# Fit a normal distribution to the data
mu_virus_nm, std_virus_nm = norm.fit(virus_nm_result7)
mu_virus_pt, std_virus_pt = norm.fit(virus_pt_result7)

mu_bacteria_nm, std_bacteria_nm = norm.fit(bacteria_nm_result7)
mu_bacteria_pt, std_bacteria_pt = norm.fit(bacteria_pt_result7)

mu_archaea_nm, std_archaea_nm = norm.fit(archaea_nm_result7)
mu_archaea_pt, std_archaea_pt = norm.fit(archaea_pt_result7)

mu_fungi_nm, std_fungi_nm = norm.fit(fungi_nm_result7)
mu_fungi_pt, std_fungi_pt = norm.fit(fungi_pt_result7)

mu_plant_nm, std_plant_nm = norm.fit(plant_nm_result7)
mu_plant_pt, std_plant_pt = norm.fit(plant_pt_result7)

mu_protozoa_nm, std_protozoa_nm = norm.fit(protozoa_nm_result7)
mu_protozoa_pt, std_protozoa_pt = norm.fit(protozoa_pt_result7)

mu_mito_nm, std_mito_nm = norm.fit(mito_nm_result7)
mu_mito_pt, std_mito_pt = norm.fit(mito_pt_result7)

mu_plastid_nm, std_plastid_nm = norm.fit(plastid_nm_result7)
mu_plastid_pt, std_plastid_pt = norm.fit(plastid_pt_result7)

# # Plot the histogram.
# plt.hist(virus_nm_result7, bins=25, density=True, alpha=0.6, color='g')

# # Plot the PDF.
# xmin, xmax = plt.xlim()
# x = np.linspace(xmin, xmax, 100)
# p_virus_nm = norm.pdf(x, mu_virus_nm, std_virus_nm)

# plt.plot(x, p_virus_nm, 'k', linewidth=2)
# title = "Viral Nucleotides: \u03BC=%.2f,  \u03C3=%.2f" % (mu_virus_nm, std_virus_nm)
# plt.title(title)
# plt.ylabel("Frequency (%)")
# plt.xlabel("Normalized Dissimilarity Rate")

# plt.show()

print(NormalDist(mu_virus_nm,std_virus_nm).pdf(0.91))
print(NormalDist(mu_virus_pt,std_virus_pt).pdf(0.91))