'''
  Usage: python3 Tests/accuracy_CV.py CV_NBG_With_LEN Accuracy
         python3 Tests/accuracy_CV.py KNN/CrossValidation Accuracy
         python3 Tests/accuracy_CV.py CV_NBG_With_LEN F1Score
         python3 Tests/accuracy_CV.py KNN/CrossValidation F1Score
'''

from sklearn.metrics import f1_score
from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
import sys

# Accuracy
Hits = 0;
TotalSamples = 0;

# F1_score
Test_Domains = []
Train_Domains = []

folder = sys.argv[1]
Mode = str(sys.argv[2])

f = open("Tests/Predictions/"+folder+"/Prediction_Viral_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Viral")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Virus' or prediction.rstrip('\n') == 'Viral':
    Hits = Hits + 1;
    
f = open("Tests/Predictions/"+folder+"/Prediction_Bacteria_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Bacteria")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Bacteria':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Archaea_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Archaea")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Archaea':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Fungi_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Fungi")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Fungi':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Plant_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Plant")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Plant':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Protozoa_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Protozoa")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Protozoa':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Mitochondrial_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Mitochondrial")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Mitochondrial':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Plastid_CV.txt", "r")
for prediction in f:
  Test_Domains.append("Plastid")
  Train_Domains.append(prediction.rstrip('\n'))
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Plastid':
    Hits = Hits + 1;

if Mode == "Accuracy":
  percentage = float(Hits) / float(TotalSamples) * 100
  print("Accuracy of : %.3f%%" % (percentage))

elif Mode == "F1Score":
  f1score=f1_score(Test_Domains, Train_Domains, average='weighted')
  print("F1 score of : %.3f%%" % (f1score * 100))

else:
  print("Mode not found. Please insert a valid mode.")