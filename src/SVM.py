
# Packages for analysis
import pandas as pd
import numpy as np
from sklearn import svm

# Packages for visuals
import matplotlib.pyplot as plt
import seaborn as sns; sns.set(font_scale=1.2)

from mpl_toolkits.mplot3d import Axes3D

samples = pd.read_csv('Analysis/SVM/Domains.csv')

samplesDomains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]

# Create a function to guess when a recipe is a Viral or a cupcake
def predictDomain(DNA, AA, GC):
    if(model.predict([[DNA, AA, GC]]))==0:
        print('You\'re looking at a Viral domain!')
    else:
        print('You\'re looking at other domain!')

# Specify inputs for the model
predictors = samples[['DNA','AA','GC']].to_numpy()
domainLabel = np.where(samples['Domain']=='Viral', 0, 1)

# Fit the SVM model
model = svm.SVC(kernel='linear', decision_function_shape='ovo')
model.fit(predictors, domainLabel)

# Get the separating hyperplane (2D)
w = model.coef_[0]
a = -w[0] / w[1]
xx = np.linspace(0, 1.5)
yy = a * xx - (model.intercept_[0]) / w[1]

# Plot the parallels to the separating hyperplane that pass through the support vectors (2D)
b = model.support_vectors_[0]
yy_down = a * xx + (b[1] - a * b[0])
b = model.support_vectors_[-1]
yy_up = a * xx + (b[1] - a * b[0])

# Get the separating hyperplane (3D)
tmp = np.linspace(0,1,10)
x_hyperplane,y_hyperplane = np.meshgrid(tmp, tmp)
z_hyperplane = lambda x,y: (-model.intercept_[0]-model.coef_[0][0]*x -model.coef_[0][1]*y) / model.coef_[0][2]

# Look at the margins and support vectors - (2D)
sns.lmplot(x='DNA', y='AA', data=samples, hue='Domain', palette='Set1', fit_reg=False, scatter_kws={"s": 70})
plt.plot(xx, yy, linewidth=2, color='black')
plt.plot(xx, yy_down, 'k--')
plt.plot(xx, yy_up, 'k--')
plt.scatter(model.support_vectors_[:, 0], model.support_vectors_[:, 1], s=80, facecolors='none');
plt.plot(0.828, 0.87, 'yo', markersize='9');


# axes instance
fig3d = plt.figure(figsize=(6,6))
ax = Axes3D(fig3d, auto_add_to_figure=False)
fig3d.add_axes(ax)

# get colormap
colors = {'Viral':'red', 'Bacteria':'green', 'Archaea':'blue', 'Fungi':'yellow', 'Plant':'purple', 'Protozoa':'pink', 'Mitochondrial':'grey', 'Plastid':'brown'}


# plot samples
sc = ax.scatter(samples['DNA'], samples['AA'], samples['GC'], s=40, c=samples['Domain'].map(colors), marker='o', alpha=1)
ax.set_xlabel("DNA Compression")
ax.set_ylabel("Protein Compression")
ax.set_zlabel("GC-Content")

# Plot the hyperplane
ax.plot_surface(x_hyperplane, y_hyperplane, z_hyperplane(x_hyperplane, y_hyperplane))


# Predict Domain
predictDomain(0.828, 0.87, 0.20996284)

# legend
for domain in samplesDomains:
    plt.scatter([], [], s=40, label=str(domain))
plt.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Domains: ', loc="lower center", borderaxespad=-10, ncol=4)

# save
plt.savefig("Analysis/SVM/domains_scatter", bbox_inches='tight', dpi=300)
plt.show()