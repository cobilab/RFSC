#!/usr/bin/env python

'''
  Usage: python3 Mutations/classification.py File Accuracy
         python3 Mutations/classification.py File F1Score
         python3 Mutations/classification.py File
'''
import warnings
from xgboost import XGBClassifier 
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import classification_report,confusion_matrix
from xgboost import XGBClassifier
from sklearn.metrics import f1_score 
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.svm import SVC
import csv
import sys
import numpy as np
import os
from sklearn.metrics import classification_report


training_path="../Analysis/KNN/"
test_path="./features/"
mutationLevels=["0","1","2","4","6","8","10"]

def warn(*args, **kwargs):
    pass
warnings.warn = warn

def ReadTrainData(test_list):
    #do not change this
    domains = {
        "Viral":0, 
        "Bacteria":1,
        "Archaea":2, 
        "Fungi":3,
        "Plant":4, 
        "Protozoa":5,
        "Mitochondrial":6, 
        "Plastid":7
    }
    
    X_train, y_train = [], []
    with open(os.path.join(training_path, "Train.csv"), 'r') as file:
        samples = csv.reader(file)
        next(samples)
        for row in samples:
            values=[float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])]
            label=domains[row[0]]
            line=[label]+values
            if line not in test_list:
                X_train.append(values)
                y_train.append(label)
    with open(os.path.join(training_path, "Test.csv"), 'r') as file:
        samples = csv.reader(file)
        next(samples)
        for row in samples:
            values=[float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])]
            label=domains[row[0]]
            line=[label]+values
            if line not in test_list:
                X_train.append(values)
                y_train.append(label)
    
    return np.array(X_train).astype('float32'), np.array(y_train).astype('int32')


def ReadTestData(filename):
    domains = {
        "viral":0, 
        "bacteria":1,
        "archaea":2, 
        "fungi":3,
        "plant":4, 
        "protozoa":5,
        "mitochondrion":6, 
        "plastid":7
    }
    
    X_test, y_test = [], []
    with open(os.path.join(test_path, filename), 'r') as file:
        samples = csv.reader(file)
        next(samples)
        for row in samples:
            X_test.append([float(row[1]),float(row[2]),float(row[3]),float(row[4]),float(row[5])])
            y_test.append(domains[row[0]])
    
    test_list=[[l]+a for l,a in zip(y_test,X_test)]  
    return np.array(X_test).astype('float32'), np.array(y_test).astype('int32'), test_list


def Classify():
    filename = "mutation_level_{}.csv".format(str(sys.argv[1]))
    if len(sys.argv) == 2:
        Mode="ALL"
    else:
        Mode = str(sys.argv[2])

    domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrion", "Plastid"]

    X_test, y_test, test_list = ReadTestData(filename)
    X_train, y_train = ReadTrainData(test_list)

    model = XGBClassifier(max_depth=12, learning_rate=0.89, n_estimators=500, eval_metric='mlogloss')
    model.fit(X_train, y_train)
    XGB_y_pred = model.predict(X_test)
    XGB_predictions = [round(value) for value in XGB_y_pred]

    GNB_model = GaussianNB()
    KNN_model = KNeighborsClassifier(n_neighbors=len(domains))
    
    GNB_model.fit(X_train, y_train)
    KNN_model.fit(X_train, y_train)
    GNB_prediction = GNB_model.predict(X_test)
    KNN_prediction = KNN_model.predict(X_test)

    if Mode == "Accuracy":
        accuracy_GNB = accuracy_score(y_test, GNB_prediction)
        accuracy_KNN = accuracy_score(y_test, KNN_prediction)
        accuracy_XGB = accuracy_score(y_test, XGB_predictions)
        print("Accuracy of GNB: %.2f%%" % (accuracy_GNB * 100.0))
        print("Accuracy of KNN: %.2f%%" % (accuracy_KNN * 100.0))
        print("Accuracy of XGB: %.2f%%" % (accuracy_XGB * 100.0))

    elif Mode == "F1Score":
        f1score_GNB=f1_score(y_test, GNB_prediction, average='weighted')
        f1score_KNN=f1_score(y_test, KNN_prediction, average='weighted')
        f1score_XGB=f1_score(y_test, XGB_predictions, average='weighted')
        print("F1 score of GNB: %.2f%%" % (f1score_GNB * 100.0))
        print("F1 score of KNN: %.2f%%" % (f1score_KNN * 100.0))
        print("F1 score of XGB: %.2f%%" % (f1score_XGB * 100.0))
    
    elif Mode == "ALL":
        print("GNB")
        print(classification_report(y_test, GNB_prediction, target_names=domains, digits=4))
        print("KNN")
        print(classification_report(y_test, KNN_prediction, target_names=domains, digits=4))
        print("XGB")
        print(classification_report(y_test, XGB_predictions, target_names=domains, digits=4))


    else:
        print("Mode not found. Please insert a valid mode. E.g.: F1Score or Accuracy")


def get_accuracy(predictions):
    accuracy = accuracy_score(y_test, predictions)
    print("Accuracy of : %.2f%%" % (accuracy * 100.0))

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of Mutations/! There are relative paths defined in this code that need to be respected!")
    elif len(sys.argv) < 2:
        print("ERROR: Missing parameter: Mutation file with features")
    elif sys.argv[1] not in mutationLevels:
        print("ERROR: Mutation level not recognized!")
    else:
        warnings.filterwarnings(action='ignore', category=DeprecationWarning)
        Classify()