import os
import sys
from os import listdir, path, makedirs
from os.path import isfile, join
sys.path.append(os.path.abspath('../src'))
from find_greatest_length import maxvaluefound_length_DNA, maxvaluefound_length_AA
from pprint import pprint
from collections import Counter

root = "./features/"
mutName = "mutation_level_"
percentages = [0,1,2,4,6,8,10]
tmpDir = os.path.join(root, "tmp/")
sequencesDir = "./original_sequences/"


def main():
    _initialize()
    for level in percentages:
        csvContent = getSequences(level)
        writeCSV(level, csvContent)
        os.system(f"rm -r {tmpDir}")

def getSequences(mutationLevel):
    content = [["Domain","DNA","AA","GC","L_DNA","L_AA"]]
    for domain in listdir(sequencesDir):
        print(f"Computing {domain}...")
        tmpPath = join(tmpDir,domain)
        if not path.exists(tmpPath):
            makedirs(tmpPath)

        for fileName in listdir(join(sequencesDir, domain)):
            name = fileName.replace(".fna.gz", "")
            csvEntry = {
                "Domain": domain,
                "DNA": None,
                "AA": None,
                "GC": None,
                "L_DNA": None,
                "L_AA": None
            }

            print(f"\tFile {fileName}")
            filePath = join(sequencesDir, domain, fileName)
            if isfile(filePath):
                if mutationLevel>0:
                    mutate(mutationLevel, tmpPath, filePath, name)
                os.system(f"../ORFs/orfm/orfm {filePath} > {tmpPath}/PROTEIN.fasta")
                os.system(f'cat {tmpPath}/PROTEIN.fasta | grep -v ">" | tr -d -c "ACDEFGHIKLMNPQRSTVWY" > {tmpPath}/{name}.fna')
                os.system(f"rm {tmpPath}/PROTEIN.fasta")

                with open(join(tmpPath, f"{name}.fna"), 'r') as fp:
                    csvEntry["L_AA"] = str(float("{:.8f}".format(len(fp.read())/maxvaluefound_length_AA)))

                os.system(f'zcat {filePath} | grep -v ">" | tr -d -c "ACGT" > {tmpPath}/GENOME_FILE')
                with open(join(tmpPath, "GENOME_FILE"), 'r') as fp:
                    c = fp.read()
                    r = Counter(c)
                    csvEntry["GC"] = str(float("{:.8f}".format((r["G"] + r["C"])/(r["G"] + r["C"] + r["A"] + r["T"]))))
                    csvEntry["L_DNA"] = str(float("{:.8f}".format(len(c)/maxvaluefound_length_DNA)))

                os.system(f'zcat {filePath} | grep -v ">" | tr -d -c "ACGT" > {tmpPath}/TO_COMPRESS')
                os.system(f"GeCo3 -v -l 3 {tmpPath}/TO_COMPRESS | sed '1,6d' | sed '2d' > {tmpPath}/RESULTS_GECO")
                
                with open(join(tmpPath, "RESULTS_GECO"), 'r') as fp:
                    try:
                        csvEntry["DNA"] = str(float(fp.read().split(" ")[7])/2)
                    except:
                        csvEntry["DNA"] = ""
                    
                os.system(f"AC -v -l 3 {tmpPath}/{name}.fna | sed '3,4d' | sed 'N;s/\\n/ /' > {tmpPath}/RESULTS_AC")
                with open(join(tmpPath, "RESULTS_AC"), 'r') as fp:
                    try:
                        csvEntry["AA"] = str(float(fp.read().split(" ")[15]))
                    except:
                        csvEntry["AA"] = ""
                os.system(f"rm {tmpPath}/*")

                content.append(csvEntry.values())
    return content

def mutate(mutationLevel, tmpPath, filePath, name):
    os.system(f"cp {filePath} {tmpPath}/{name}.fna.gz")
    os.system(f"gunzip {tmpPath}/{name}.fna.gz")
    os.system(f"gto_fasta_mutate -e {mutationLevel} -a < {tmpPath}/{name}.fna > {tmpPath}/Temporary.fna")
    os.system(f"rm {tmpPath}/{name}.fna")
    os.system(f"mv {tmpPath}/Temporary.fna {tmpPath}/{name}.fna")
    os.system(f"gzip {tmpPath}/{name}.fna")
    

def writeCSV(mutationLevel, content):
    fileName = mutName+str(mutationLevel)+".csv"
    f = open(join(root, fileName), 'w')
    for line in content:
        f.write(",".join(line))
        f.write("\n")
    f.close()

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