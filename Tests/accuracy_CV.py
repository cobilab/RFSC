'''
  Usage: python3 Tests/accuracy_CV.py CV_NBG_With_LEN
         python3 Tests/accuracy_CV.py KNN/CrossValidation
'''

from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
import sys

Hits = 0;
TotalSamples = 0;

folder = sys.argv[1]

f = open("Tests/Predictions/"+folder+"/Prediction_Viral_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Virus' or prediction.rstrip('\n') == 'Viral':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Bacteria_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Bacteria':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Archaea_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Archaea':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Fungi_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Fungi':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Plant_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Plant':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Protozoa_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Protozoa':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Mitochondrial_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Mitochondrial':
    Hits = Hits + 1;

f = open("Tests/Predictions/"+folder+"/Prediction_Plastid_CV.txt", "r")
for prediction in f:
  TotalSamples = TotalSamples + 1;
  if prediction.rstrip('\n') == 'Plastid':
    Hits = Hits + 1;

percentage = float(Hits) / float(TotalSamples) * 100

print("Percentage of success: " + str("{:.2f}".format(percentage)) + "%")