import os
import sys
from os import listdir, path, makedirs
from os.path import isfile, join
sys.path.append(os.path.abspath('../src'))
from find_greatest_length import maxvaluefound_length_DNA, maxvaluefound_length_AA
from pprint import pprint

root = "./features/"
mutName = "mutation_level_"
percentages = [0,1,2,4,6,8,10]
tmpDir = os.path.join(root, "tmp/")
#domains = ["viral", "bacteria", "archaea", "fungi", "plant", "protozoa", "mitochondrion", "plastid"]
sequencesDir = "./original_sequences/"


def main():
    _initialize()
    sequences = getSequences()

def getSequences():
    for domain in listdir(sequencesDir):
        print(f"Computing {domain}...")
        tmpPath = join("tmp",domain,"GeneratedORFs")

        print(tmpPath)
        if not path.exists(tmpPath):
            makedirs(tmpPath)

        for fileName in listdir(join(sequencesDir, domain)):
            name = fileName.replace(".fna.gz", "")
            print(f"\tFile {fileName}")
            filePath = join(sequencesDir, domain, fileName)
            if isfile(filePath):
                os.system(f"../ORFs/orfm/orfm {filePath} > {tmpPath}/PROTEIN.fasta")
                print(name, fileName)
                os.system(f'cat {tmpPath}/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > {tmpPath}/{name}.fna')
            break# remove this
        break# remove this


def _initialize():
    if not path.exists(root):
        makedirs(root)
    if not path.exists(tmpDir):
        makedirs(tmpDir)

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of Mutations/! There are relative paths defined in this code that need to be respected!")
    else:
        main()