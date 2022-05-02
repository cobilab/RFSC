import sys
from os import listdir, path, makedirs
from os.path import isfile, join
from pprint import pprint
import os
import warnings

dbLocation="./krakenDB"
synth_sequences="../synth_merged_sequences/"

def compute():
    right = {}
    wrong = {}
    for domain in listdir(synth_sequences):
        right[domain] = 0
        wrong[domain] = 0
        for file in listdir(join(synth_sequences, domain)):
            fileName=file.split(".fna")[0]
            reportName=join("./report/", fileName)
            sequence=join(synth_sequences, domain, file)
            os.system(f"kraken2 --db {dbLocation} --report {reportName} --report-minimizer-data --output output.txt {sequence}")
            with open(reportName, 'r') as fp:
                content = fp.read().rstrip().lower()
                if domain.lower() in content:
                    right[domain] += 1
                else:
                    wrong[domain] += 1
        pprint(right)
        pprint(wrong)
    for domain in listdir(synth_sequences):
        rights = right[domain]
        total = right[domain] + wrong[domain]
        acc = rights/total
        print(f"Accuracy for {domain}: {acc}")

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of Kraken2/! There are relative paths defined in this code that need to be respected!")
    else:
        warnings.filterwarnings(action='ignore', category=DeprecationWarning)
        compute()
