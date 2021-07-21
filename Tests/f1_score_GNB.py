'''
  Usage: python3 Tests/fi_score_GNB.py
'''

from matplotlib import font_manager
import numpy as np
from scipy.stats import norm
import matplotlib.pyplot as plt
import csv
import sys

TP = [0,0,0,0,0,0,0,0]
FP = [0,0,0,0,0,0,0,0]
FN = [0,0,0,0,0,0,0,0]
TN = [0,0,0,0,0,0,0,0]

f = open("Tests/Predictions/KNN/TestDatabase.txt", "r")
for prediction in f:
    if prediction.rstrip('\n').split(",")[0] == 'Viral':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[0] = TP[0] + 1
        else:
            FN[0] = FN[0] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Viral' and prediction.rstrip('\n').split(",")[1] == 'Viral':
        FP[0] = FP[0] + 1
  
    if prediction.rstrip('\n').split(",")[0] == 'Bacteria':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[1] = TP[1] + 1
        else:
            FN[1] = FN[1] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Bacteria' and prediction.rstrip('\n').split(",")[1] == 'Bacteria':
        FP[1] = FP[1] + 1

    if prediction.rstrip('\n').split(",")[0] == 'Archaea':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[2] = TP[2] + 1
        else:
            FN[2] = FN[2] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Archaea' and prediction.rstrip('\n').split(",")[1] == 'Archaea':
        FP[2] = FP[2] + 1

    if prediction.rstrip('\n').split(",")[0] == 'Fungi':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[3] = TP[3] + 1
        else:
            FN[3] = FN[3] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Fungi' and prediction.rstrip('\n').split(",")[1] == 'Fungi':
        FP[3] = FP[3] + 1

    if prediction.rstrip('\n').split(",")[0] == 'Plant':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[4] = TP[4] + 1
        else:
            FN[4] = FN[4] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Plant' and prediction.rstrip('\n').split(",")[1] == 'Plant':
        FP[4] = FP[4] + 1

    if prediction.rstrip('\n').split(",")[0] == 'Protozoa':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[5] = TP[5] + 1
        else:
            FN[5] = FN[5] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Protozoa' and prediction.rstrip('\n').split(",")[1] == 'Protozoa':
        FP[5] = FP[5] + 1

    if prediction.rstrip('\n').split(",")[0] == 'Mitochondrial':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[6] = TP[6] + 1
        else:
            FN[6] = FN[6] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Mitochondrial' and prediction.rstrip('\n').split(",")[1] == 'Mitochondrial':
        FP[6] = FP[6] + 1

    if prediction.rstrip('\n').split(",")[0] == 'Plastid':
        if prediction.rstrip('\n').split(",")[0] == prediction.rstrip('\n').split(",")[1]:
            TP[7] = TP[7] + 1
        else:
            FN[7] = FN[7] + 1
    elif prediction.rstrip('\n').split(",")[0] != 'Plastid' and prediction.rstrip('\n').split(",")[1] == 'Plastid':
        FP[7] = FP[7] + 1


Precision = sum(TP) / (sum(TP) + sum(FP))
Recall = sum(TP) / (sum(TP) + sum(FN))
F_Score = 2 * (( Precision * Recall ) / ( Precision + Recall ))

print(F_Score)