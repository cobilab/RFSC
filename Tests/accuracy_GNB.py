'''
  Usage: python3 Tests/accuracy_GNB.py GNB/CrossValidation 0.8 Accuracy
         python3 Tests/accuracy_GNB.py GNB/CrossValidation 0.8 F1Score
'''

from sklearn.metrics import f1_score
from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import numpy as np
import os, os.path
import csv
import sys

# Accuracy
Hits = 0;
TotalSamples = 0;

# F1_score
Test_Domains = []
Train_Domains = []

folder = sys.argv[1]
trainDatabase = float(sys.argv[2])  # Percentage of the (current) domain database used for training the model (80%)
Mode = str(sys.argv[3])

dir_path = os.path.dirname(os.path.realpath(__file__))  # Current directory path

virusFile = dir_path + '/../Analysis/GeCo/Virus/geco3_Viral.csv'
bacteriaFile = dir_path + '/../Analysis/GeCo/Bacteria/geco3_Bacteria.csv'
archaeaFile = dir_path + '/../Analysis/GeCo/Archaea/geco3_Archaea.csv'
fungiFile = dir_path + '/../Analysis/GeCo/Fungi/geco3_Fungi.csv'
plantFile = dir_path + '/../Analysis/GeCo/Plant/geco3_Plant.csv'
protozoaFile = dir_path + '/../Analysis/GeCo/Protozoa/geco3_Protozoa.csv'
mitoFile = dir_path + '/../Analysis/GeCo/Mitochondrial/geco3_Mitochondrial.csv'
plastidFile = dir_path + '/../Analysis/GeCo/Plastid/geco3_Plastid.csv'

f = open("Tests/Predictions/"+folder+"/Prediction_Viral.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Virus' or prediction.rstrip('\n') == 'Viral':
    Hits = Hits + 1;

numVirusFiles = int(len(open(virusFile).readlines(  )) * trainDatabase)
Train_Domains = np.full(numVirusFiles,"Viral")

    
f = open("Tests/Predictions/"+folder+"/Prediction_Bacteria.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Bacteria':
    Hits = Hits + 1;

numBacteriaFiles = int(len(open(bacteriaFile).readlines(  )) * trainDatabase)
Train_Bacteria = np.full(numBacteriaFiles,"Bacteria")
Train_Domains = np.concatenate((Train_Domains, Train_Bacteria), axis=None)

f = open("Tests/Predictions/"+folder+"/Prediction_Archaea.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Archaea':
    Hits = Hits + 1;

numArchaeaFiles = int(len(open(archaeaFile).readlines(  )) * trainDatabase)
Train_Archaea = np.full(numArchaeaFiles,"Archaea")
Train_Domains = np.concatenate((Train_Domains, Train_Archaea), axis=None)

f = open("Tests/Predictions/"+folder+"/Prediction_Fungi.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Fungi':
    Hits = Hits + 1;

numFungiFiles = int(len(open(fungiFile).readlines(  )) * trainDatabase)
Train_Fungi = np.full(numFungiFiles,"Fungi")
Train_Domains = np.concatenate((Train_Domains, Train_Fungi), axis=None)

f = open("Tests/Predictions/"+folder+"/Prediction_Plant.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Plant':
    Hits = Hits + 1;

numPlantFiles = int(len(open(plantFile).readlines(  )) * trainDatabase)
Train_Plant = np.full(numPlantFiles,"Plant")
Train_Domains = np.concatenate((Train_Domains, Train_Plant), axis=None)

f = open("Tests/Predictions/"+folder+"/Prediction_Protozoa.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Protozoa':
    Hits = Hits + 1;

numProtozoaFiles = int(len(open(protozoaFile).readlines(  )) * trainDatabase)
Train_Protozoa = np.full(numProtozoaFiles,"Protozoa")
Train_Domains = np.concatenate((Train_Domains, Train_Protozoa), axis=None)

f = open("Tests/Predictions/"+folder+"/Prediction_Mitochondrial.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Mitochondrial':
    Hits = Hits + 1;

numMitoFiles = int(len(open(mitoFile).readlines(  )) * trainDatabase)
Train_Mito = np.full(numMitoFiles,"Mitochondrial")
Train_Domains = np.concatenate((Train_Domains, Train_Mito), axis=None)

f = open("Tests/Predictions/"+folder+"/Prediction_Plastid.txt", "r")
for prediction in f:
  Test_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Plastid':
    Hits = Hits + 1;

numPlastidFiles = int(len(open(plastidFile).readlines(  )) * trainDatabase)
Train_Plastid = np.full(numPlastidFiles,"Plastid")
Train_Domains = np.concatenate((Train_Domains, Train_Plastid), axis=None)

if Mode == "Accuracy":
  percentage = float(Hits) / float(TotalSamples) * 100
  print("Accuracy of : %.3f%%" % (percentage))

elif Mode == "F1Score":
  f1score=f1_score(Test_Domains, Train_Domains.tolist(), average='weighted')
  print("F1 score of : %.3f%%" % (f1score * 100))

else:
  print("Mode not found. Please insert a valid mode.")