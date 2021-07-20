'''
  Usage: python3 Tests/accuracy_KNN.py
'''

from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
import sys

Hits = [0,0,0,0,0,0,0,0];
Miss = [0,0,0,0,0,0,0,0];
TotalPredictions = 0;

f = open("Tests/Predictions/KNN/TestDatabase.txt", "r")
for prediction in f:
  #print(prediction.rstrip('\n').split(",")[0])
  if prediction.rstrip('\n').split(",")[0] == 'Viral':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[0] = Hits[0] + 1;
    else:
      Miss[0] = Miss[0] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Bacteria':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[1] = Hits[1] + 1;
    else:
      Miss[1] = Miss[1] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Archaea':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[2] = Hits[2] + 1;
    else:
      Miss[2] = Miss[2] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Fungi':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[3] = Hits[3] + 1;
    else:
      Miss[3] = Miss[3] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Plant':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[4] = Hits[4] + 1;
    else:
      Miss[4] = Miss[4] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Protozoa':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[5] = Hits[5] + 1;
    else:
      Miss[5] = Miss[5] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Mitochondrial':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[6] = Hits[6] + 1;
    else:
      Miss[6] = Miss[6] + 1;
  elif prediction.rstrip('\n').split(",")[0] == 'Plastid':
    if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
      Hits[7] = Hits[7] + 1;
    else:
      Miss[7] = Miss[7] + 1;

TotalPredictions = sum(Hits) + sum (Miss)

print(Hits)
print(Miss)
print(TotalPredictions)

print(sum(Hits) / TotalPredictions * 100)