#!/usr/bin/env python
import os, os.path, sys
from os.path import join
from os import listdir, path, makedirs
import warnings

synth_sequences="./Kraken2/report"
domains =["viral","bacteria","archaea", "fungi","plant", "protozoa","mitochondrion", "plastid"]
def compute():
    right = {}
    wrong = {}
    for domain in domains:
        right[domain] = 0
        wrong[domain] = 0
    for file in listdir(synth_sequences):
        # fileName=file.split(".fna")[0]
        fileName = ''.join([i for i in file if not i.isdigit()])
        reportName=join(synth_sequences, file)
        with open(reportName, 'r') as fp:
            content = fp.read().rstrip().lower()
            if fileName.lower() in content:
                right[fileName] += 1
            else:
                wrong[fileName] += 1
    accuracy={}
    for domain in domains:
        accuracy[domain]=(right[domain]/(right[domain]+wrong[domain]))*100
    accuracy["average"]=(sum(list(accuracy.values()))/len(list(accuracy.values())))

    print("Accuracy:")
    for i in accuracy:
        print(i.title(),":",accuracy[i])

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of Kraken2/! There are relative paths defined in this code that need to be respected!")
    else:
        warnings.filterwarnings(action='ignore', category=DeprecationWarning)
        compute()