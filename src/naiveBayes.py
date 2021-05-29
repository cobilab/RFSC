import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
from statistics import NormalDist

# ( https://www.youtube.com/watch?v=rzFX5NWojp0 ) 

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

# Fit a normal distribution to the data:
mu, std = norm.fit(virus_nm_result7)

# Plot the histogram.
plt.hist(virus_nm_result7, bins=25, density=True, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)

plt.show()

#print(x)
#print("-----")
#print(p)
#print("-----")
print(NormalDist(mu,std).pdf(0.91))