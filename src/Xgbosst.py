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

def warn(*args, **kwargs):
    pass
warnings.warn = warn

def predict_XGBClassifier():
    domains = ["Viral", "Bacteria", "Archaea", "Fungi", "Plant", "Protozoa", "Mitochondrial", "Plastid"]

    Mode = str(sys.argv[1])				# Test | CV

    crossValidation_Min_Lim = float(sys.argv[2])    # Inferior limit of the block used for testing
    crossValidation_Max_Lim = float(sys.argv[3])    # Superior limit of the block used for testing

    dna_input = float(sys.argv[4])   	# DNA Compression Value
    aa_input = float(sys.argv[5])    	# AA Compression Value
    gc_content = float(sys.argv[6])    	# GC-Content Value
    len_dna_input = float(sys.argv[7])  # Length of the DNA sequence
    len_aa_input = float(sys.argv[8])   # Length of the AA sequence

    CV_Domain = str(sys.argv[9])		# NULL | Viral

    X_train, y_train = [], []
    X_test = []

    X_test.append([dna_input, aa_input, gc_content, len_dna_input, len_aa_input])

    if Mode == "Test":
        with open('Analysis/Xgbosst/Train.csv', 'r') as file:
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

    elif Mode == "CV":
        numberSamples = [0,0,0,0,0,0,0,0,0]
        with open('Analysis/Xgbosst/Train.csv', 'r') as file:
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

        with open('Analysis/Xgbosst/Train.csv', 'r') as file:
            samples = csv.reader(file)
            next(samples)
            for row in samples:
                if sampleCounter >= numFiles_ForTestMin and sampleCounter <= numFiles_ForTestMax:
                    sampleCounter = sampleCounter + 1
                else:
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

    else:
    	print("Mode not found. Please specify the correct mode of use: Test or CV")
    
    X_train = np.array(X_train).astype('float32')
    y_train = np.array(y_train).astype('int32') 
    X_test = np.array(X_test).astype('float32')

    model = XGBClassifier(max_depth=12, learning_rate=0.89, n_estimators=500, eval_metric='mlogloss')
    model.fit(X_train, y_train)
    prediction = model.predict(X_test)
    print(domains[prediction[0]])

if __name__ == "__main__":
    predict_XGBClassifier()