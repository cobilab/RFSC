import sys
from os import listdir, path, makedirs
from os.path import isfile, join

dbLocation="./krakenDB"
synth_sequences="../synth_merged_sequences/"

def compute():
    for domain in listdir(synth_sequences):
        for file in listdir(join(synth_sequences, domain)):
            reportName=file.split(".fna")[0]
            reportName=join("./report/", reportName)
            sequence=join(synth_sequences, domain, file)
            os.system(f"kraken2 --db {dbLocation} --report {reportName} --report-minimizer-data --output output.txt {sequence}")

def classify():
    print("to do")

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of Kraken2/! There are relative paths defined in this code that need to be respected!")
    else:
        warnings.filterwarnings(action='ignore', category=DeprecationWarning)
        compute()
        classify()
