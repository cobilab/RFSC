import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
from statistics import NormalDist

virus_nm_size1, virus_nm_size2, virus_nm_size3, virus_nm_size4, virus_nm_size7, virus_nm_size12, virus_nm_size14 = [], [], [], [], [], [], []
virus_nm_result1, virus_nm_result2, virus_nm_result3, virus_nm_result4, virus_nm_result7, virus_nm_result12, virus_nm_result14 = [], [], [], [], [], [], []
virus_pt_size1, virus_pt_size2, virus_pt_size3, virus_pt_size4, virus_pt_size7, virus_pt_size12, virus_pt_size14 = [], [], [], [], [], [], []
virus_pt_result1, virus_pt_result2, virus_pt_result3, virus_pt_result4, virus_pt_result7, virus_pt_result12, virus_pt_result14 = [], [], [], [], [], [], []


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

        virus_nm_result1.append("0"+row[0].split("\t")[8])
        virus_nm_result2.append("0"+row[0].split("\t")[9])
        virus_nm_result3.append("0"+row[0].split("\t")[10])
        virus_nm_result4.append("0"+row[0].split("\t")[11])
        virus_nm_result7.append("0"+row[0].split("\t")[12])
        virus_nm_result12.append("0"+row[0].split("\t")[13])
        virus_nm_result14.append("0"+row[0].split("\t")[14])

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

        virus_pt_result1.append("0"+row[0].split("\t")[8])
        virus_pt_result2.append("0"+row[0].split("\t")[9])
        virus_pt_result3.append("0"+row[0].split("\t")[10])
        virus_pt_result4.append("0"+row[0].split("\t")[11])
        virus_pt_result7.append("0"+row[0].split("\t")[12])
        virus_pt_result12.append("0"+row[0].split("\t")[13])
        virus_pt_result14.append("0"+row[0].split("\t")[14])

# Generate some data for this demonstration.
#data = norm.rvs(10.0, 2.5, size=500)
#data = [8, 12, 7, 6, 2, 6, 7, 4, 3, 2, 5]
virus_nm_result1 = np.array(virus_nm_result1[1:], dtype=float)
virus_nm_result2 = np.array(virus_nm_result2[1:], dtype=float)
virus_nm_result3 = np.array(virus_nm_result3[1:], dtype=float)
virus_nm_result4 = np.array(virus_nm_result4[1:], dtype=float)
virus_nm_result7 = np.array(virus_nm_result7[1:], dtype=float)
virus_nm_result12 = np.array(virus_nm_result12[1:], dtype=float)
virus_nm_result14 = np.array(virus_nm_result14[1:], dtype=float)

# Fit a normal distribution to the data:
mu, std = norm.fit(virus_nm_result3)

# Plot the histogram.
plt.hist(virus_nm_result3, bins=25, density=True, alpha=0.6, color='g')

# Plot the PDF.
xmin, xmax = plt.xlim()
x = np.linspace(xmin, xmax, 100)
p = norm.pdf(x, mu, std)
plt.plot(x, p, 'k', linewidth=2)
title = "Fit results: mu = %.2f,  std = %.2f" % (mu, std)
plt.title(title)

plt.show()

print(x)
print("-----")
print(p)
print("-----")
print(NormalDist(mu,std).pdf(0.9))