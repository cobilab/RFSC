'''
	Usage: python3 src/KNN.py Test 2 0 1 0.7545 0.783272 0.19727258 1.927032e-10 3.87563e-11 NULL
'''

# Packages for analysis
import pandas as pd
import numpy as np
from math import sqrt
import csv, sys
import os, os.path

Mode = str(sys.argv[1])				# Test | CV

K = int(sys.argv[2]);				# K-Neighbors

crossValidation_Min_Lim = float(sys.argv[3])    # Inferior limit of the block used for testing
crossValidation_Max_Lim = float(sys.argv[4])    # Superior limit of the block used for testing

dna_input = float(sys.argv[5])   	# DNA Compression Value
aa_input = float(sys.argv[6])    	# AA Compression Value
gc_content = float(sys.argv[7])    	# GC-Content Value
len_dna_input = float(sys.argv[8])  # Length of the DNA sequence
len_aa_input = float(sys.argv[9])   # Length of the AA sequence

CV_Domain = str(sys.argv[10])		# NULL | Viral

samplesDomains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]
train_samples = []
to_predict = [dna_input, aa_input, gc_content, len_dna_input, len_aa_input]

dir_path = os.path.dirname(os.path.realpath(__file__))  # Current directory path

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

if Mode == "Test":
	with open('Analysis/KNN/Train.csv', 'r') as file:
		samples = csv.reader(file)
		next(samples)
		for row in samples:
			train_samples.append([float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5]),row[0]])

elif Mode == "CV":
	numberSamples = [0,0,0,0,0,0,0,0,0]
	with open('Analysis/KNN/Train.csv', 'r') as file:
		samples = csv.reader(file)
		next(samples)
		for row in samples:
			if row[0] == "Viral":
				numberSamples[0] = numberSamples[0] + 1;
			elif row[0] == "Bacteria":
				numberSamples[1] = numberSamples[1] + 1;
			elif row[0] == "Archaea":
				numberSamples[2] = numberSamples[2] + 1;
			elif row[0] == "Fungi":
				numberSamples[3] = numberSamples[3] + 1;
			elif row[0] == "Plant":
				numberSamples[4] = numberSamples[4] + 1;
			elif row[0] == "Protozoa":
				numberSamples[5] = numberSamples[5] + 1;
			elif row[0] == "Mitochondrial":
				numberSamples[6] = numberSamples[6] + 1;
			elif row[0] == "Plastid":
				numberSamples[7] = numberSamples[7] + 1;

	if CV_Domain == "Viral":
		limSamples = numberSamples[0]
	elif CV_Domain == "Bacteria":
		limSamples = numberSamples[0] + numberSamples[1]
	elif CV_Domain == "Archaea":
		limSamples = numberSamples[0] + numberSamples[1] + numberSamples[2]
	elif CV_Domain == "Fungi":
		limSamples = numberSamples[0] + numberSamples[1] + numberSamples[2] + numberSamples[3]
	elif CV_Domain == "Plant":
		limSamples = numberSamples[0] + numberSamples[1] + numberSamples[2] + numberSamples[3] + numberSamples[4]
	elif CV_Domain == "Protozoa":
		limSamples = numberSamples[0] + numberSamples[1] + numberSamples[2] + numberSamples[3] + numberSamples[4] + numberSamples[5]
	elif CV_Domain == "Mitochondrial":
		limSamples = numberSamples[0] + numberSamples[1] + numberSamples[2] + numberSamples[3] + numberSamples[4] + numberSamples[5] + numberSamples[6] 
	elif CV_Domain == "Plastid":
		limSamples = numberSamples[0] + numberSamples[1] + numberSamples[2] + numberSamples[3] + numberSamples[4] + numberSamples[5] + numberSamples[6] + numberSamples[7]

	numFiles_ForTestMin = int(limSamples * crossValidation_Min_Lim)
	numFiles_ForTestMax = int(limSamples * crossValidation_Max_Lim)

	sampleCounter = 0;

	with open('Analysis/KNN/Train.csv', 'r') as file:
		samples = csv.reader(file)
		next(samples)
		for row in samples:
			if sampleCounter >= numFiles_ForTestMin and sampleCounter <= numFiles_ForTestMax:
				sampleCounter = sampleCounter + 1
			else:
				train_samples.append([float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5]),row[0]])
				sampleCounter = sampleCounter + 1

else:
	print("Mode not found. Please specify the correct mode of use: Test or CV")
      
#for row in train_samples:
#	distance = euclidean_distance(to_predict, row)
#	print(distance)

#neighbors = get_neighbors(train_samples, to_predict, K)
#for neighbor in neighbors:
#	print(neighbor)

prediction = predict_classification(train_samples, to_predict, K)
print(prediction)