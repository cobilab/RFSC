#!/usr/bin/env python

'''
  Usage: python3 Tests/accuracy_Xgbosst.py Accuracy
         python3 Tests/accuracy_Xgbosst.py F1Score
         python3 Tests/accuracy_Xgbosst.py
'''

import warnings
from xgboost import XGBClassifier 
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report,confusion_matrix
from xgboost import XGBClassifier
from sklearn.metrics import f1_score 
import csv
import sys
import numpy as np
from sklearn.metrics import classification_report

def warn(*args, **kwargs):
    pass
warnings.warn = warn

def predict_XGBClassifier():
    domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]
    
    if len(sys.argv) > 1:
        Mode = str(sys.argv[1])
    else:
        Mode="ALL"

    X_train, y_train = [], []
    X_test, y_test = [], []

    with open('Analysis/KNN/Train.csv', 'r') as file:
        samples = csv.reader(file)
        next(samples)
        for row in samples:
            X_train.append([float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])])
            if row[0] == "Viral":
                y_train.append(0)
            elif row[0] == "Bacteria":
                y_train.append(1)
            elif row[0] == "Archaea":
                y_train.append(2)
            elif row[0] == "Fungi":
                y_train.append(3)
            elif row[0] == "Plant":
                y_train.append(4)
            elif row[0] == "Protozoa":
                y_train.append(5)
            elif row[0] == "Mitochondrial":
                y_train.append(6)
            elif row[0] == "Plastid":
                y_train.append(7)

    with open('Analysis/KNN/Test.csv', 'r') as file:
        samples = csv.reader(file)
        next(samples)
        for row in samples:
            X_test.append([float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])])
            if row[0] == "Viral":
                y_test.append(0)
            elif row[0] == "Bacteria":
                y_test.append(1)
            elif row[0] == "Archaea":
                y_test.append(2)
            elif row[0] == "Fungi":
                y_test.append(3)
            elif row[0] == "Plant":
                y_test.append(4)
            elif row[0] == "Protozoa":
                y_test.append(5)
            elif row[0] == "Mitochondrial":
                y_test.append(6)
            elif row[0] == "Plastid":
                y_test.append(7)

    X_train = np.array(X_train).astype('float32')
    y_train = np.array(y_train).astype('int32') 
    X_test = np.array(X_test).astype('float32')
    y_test = np.array(y_test).astype('int32') 

    model = XGBClassifier(max_depth=12, learning_rate=0.89, n_estimators=500, eval_metric='mlogloss')
    model.fit(X_train, y_train)
    y_pred = model.predict(X_test)
    predictions = [round(value) for value in y_pred]

    if Mode == "Accuracy":
        accuracy = accuracy_score(y_test, predictions)
        print("Accuracy of : %.2f%%" % (accuracy * 100.0))

    elif Mode == "F1Score":
        f1score=f1_score(y_test, y_pred, average='weighted')
        print("F1 score of : %.2f%%" % (f1score * 100.0))
    
    elif Mode == "ALL":
        print(classification_report(y_test, y_pred, target_names=domains, digits=4))


    else:
        print("Mode not found. Please insert a valid mode. E.g.: F1Score or Accuracy")


if __name__ == "__main__":
    warnings.filterwarnings(action='ignore', category=DeprecationWarning)
    predict_XGBClassifier()