# Packages for analysis
import pandas as pd
import numpy as np
from math import sqrt
import csv, sys

K = 3;

dna_input = float(sys.argv[1])    # 
aa_input = float(sys.argv[1])    # 
gc_input = float(sys.argv[1])    # 

#samples = pd.read_csv('Analysis/KNN/Domains.csv', header=0)

samplesDomains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]

train_samples = []

# calculate the Euclidean distance between two vectors
def euclidean_distance(row1, row2):
	distance = 0.0
	for i in range(len(row1)-1):
		distance += (row1[i] - row2[i])**2
	return sqrt(distance)

# Locate the most similar neighbors
def get_neighbors(train, test_row, num_neighbors):
	distances = list()
	for train_row in train:
		dist = euclidean_distance(test_row, train_row)
		distances.append((train_row, dist))
	distances.sort(key=lambda tup: tup[1])
	neighbors = list()
	for i in range(num_neighbors):
		neighbors.append(distances[i][0])
	return neighbors

# Make a classification prediction with neighbors
def predict_classification(train, test_row, num_neighbors):
	neighbors = get_neighbors(train, test_row, num_neighbors)
	output_values = [row[-1] for row in neighbors]
	prediction = max(set(output_values), key=output_values.count)
	return prediction


with open('Analysis/KNN/Domains.csv', 'r') as file:
    samples = csv.reader(file)
    next(samples)
    for row in samples:
        #if (row[0] == "Viral"):
            #train_samples.append([float(row[1]),float(row[2]),row[0]])
        train_samples.append([float(row[1]),float(row[2]),row[0]])
            

row0 = [dna_input, aa_input, gc_input]
#for row in train_samples:
#	distance = euclidean_distance(row0, row)
#	print(distance)

#neighbors = get_neighbors(train_samples, row0, K)
#for neighbor in neighbors:
#	print(neighbor)

prediction = predict_classification(train_samples, row0, K)
print(prediction)