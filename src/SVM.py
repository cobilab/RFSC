# Packages for analysis
import pandas as pd
import numpy as np
from sklearn import svm

# Packages for visuals
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(font_scale=1.2)

from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import ListedColormap

# Allows charts to appear in the notebook
#%matplotlib inline

# Pickle package
import pickle

domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]

# Read in domains data [NC, AA, GC]
samples_data = pd.read_csv('Analysis/SVM/Domains.csv')
#samples_data.head()

# Plot NC & AA
#sns.lmplot(x='NC', y='AA', data=samples_data, hue='Domain', palette='Set1', fit_reg=False, scatter_kws={"s": 40});
#plt.show()

# Specify inputs for the model
# .as_matrix() is depriciated; use .to_numpy() instead
predictors = samples_data[['NC','AA','GC']].to_numpy()
type_label = np.where(samples_data['Domain']=='Viral', 0, 1)

# Feature names
domain_features = samples_data.columns.values[1:].tolist()

# Fit the SVM model
model = svm.SVC(kernel='linear')
model.fit(predictors, type_label)

''' 
    Visualize the results
'''

# Get the separating hyperplane
w = model.coef_[0]
a = -w[0] / w[1]
xx = np.linspace(30, 60)
yy = a * xx - (model.intercept_[0]) / w[1]



''' 
    3D scatter for NC, AA and GC
'''
# Get data
x = samples_data['NC']
y = samples_data['AA']
z = samples_data['GC']

# axes instance
fig = plt.figure(figsize=(6,6))
ax = Axes3D(fig, auto_add_to_figure=False)
fig.add_axes(ax)

# get colormap from seaborn
cmap = ListedColormap(sns.color_palette("husl", 256).as_hex())

# plot
sc = ax.scatter(x, y, z, s=40, c=x, marker='o', cmap=cmap, alpha=1)
ax.set_xlabel("DNA Compression")
ax.set_ylabel("Protein Compression")
ax.set_zlabel("GC-Content")

''' 
    Hyperline
'''
# Plot the hyperplane


# legend
for domain in domains:
    plt.scatter([], [], s=40, label=str(domain))
plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Domains: ', loc="lower center", borderaxespad=-10, ncol=4)

# save
plt.savefig("Analysis/SVM/domains_scatter", bbox_inches='tight', dpi=300)
plt.show()