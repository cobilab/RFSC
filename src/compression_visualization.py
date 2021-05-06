import numpy as np
import matplotlib.pyplot as plt
import csv

virus_size3 = []
virus_result3 = []

with open('GeCo3_Output/Virus/geco3_Viral.csv', 'r') as file:
    reader = csv.reader(file)
    for row in reader:
        virus_size3.append(row[0].split("\t")[3])
        virus_result3.append("0"+row[0].split("\t")[10])


# Create data
N = 60
virus = (0.6 + 0.6 * np.random.rand(N), np.random.rand(N))
bacteria = (0.4+0.3 * np.random.rand(N), 0.5*np.random.rand(N))
archaea = (0.3*np.random.rand(N),0.3*np.random.rand(N))
fungi = (0.3+0.1 *np.random.rand(N),0.3*np.random.rand(N))
plant = (0.3+0.2 *np.random.rand(N),0.9*np.random.rand(N))
protozoa = (0.3+0.9 *np.random.rand(N),0.8*np.random.rand(N))
mitochondrion = (0.3+0.7 *np.random.rand(N),0.1*np.random.rand(N))
plastid = (0.3+0.8 *np.random.rand(N),0.4*np.random.rand(N))

data = (virus, bacteria, archaea, fungi, plant, protozoa, mitochondrion, plastid)
colors = ("blue", "green", "red", "purple", "orange", "grey", "pink", "yellow")
groups = ("Viruses", "Bacteria", "Archaea", "Fungi", "Plants", "Protozoa", "Mitochondrion", "Plastid")

# Create plot
fig = plt.figure()
ax = fig.add_subplot(1, 1, 1, facecolor="1.0")

for data, color, group in zip(data, colors, groups):
    x, y = data
    ax.scatter(x, y, alpha=0.8, c=color, edgecolors='none', s=30, label=group)

plt.title('DNA compression plot by domain')
plt.legend(loc=2)
plt.show()