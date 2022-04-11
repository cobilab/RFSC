#!/usr/bin/env python
'''
  Usage: python3 SyntheticSequences/computeSyntheticFeatures.py
'''


import os
import sys
from os import listdir, path, makedirs
from os.path import isfile, join
sys.path.append(os.path.abspath('../src'))
from find_greatest_length import maxvaluefound_length_DNA, maxvaluefound_length_AA
from pprint import pprint
from collections import Counter
import itertools


sequencesDir = "./original_sequences/"
mergedPath="./synth_merged_sequences"
root = "./features/"
sodPath="../SoD/"

synthName = "synthetic_features"
percentages = [0,1,2,4,6,8,10]
tmpDir = os.path.join(root, "tmp/")
tmpSeqPath = join(".","tmpSq")
domains = ["viral", "bacteria", "archaea", "fungi", "plant", "protozoa", "mitochondrion", "plastid"]


def main():
    _initialize()
    createSyntheticSequences()
    csvContent = getSequences()


def createSyntheticSequences():
    
    for x in os.walk(sequencesDir):
        if x[1]:
            for y in x[1]:
                if not path.exists(join(mergedPath,y)):
                    makedirs(join(mergedPath,y))
        if x[2]:
            merge_list=[join(x[0],z) for z in x[2]]
            toMerge=list(itertools.combinations(merge_list,2))
            mergeAsSynthetic(toMerge)

def mergeAsSynthetic(toMergeList):
    counter=0
    for pairToMerge in toMergeList:
        get_domain_string=pairToMerge[0].split("/")
        domain=get_domain_string[2]
        filename=domain+str(counter)+".fna"
        os.system(f'zcat {pairToMerge[0]} | grep -v ">" | tr -d -c "ACGT" > {tmpSeqPath}/seqA')
        os.system(f'zcat {pairToMerge[1]} | grep -v ">" | tr -d -c "ACGT" > {tmpSeqPath}/seqB')
        with open(join(tmpSeqPath, "seqA"), 'r') as fp:
            len_A=float("{:.8f}".format(len(fp.read())))
        with open(join(tmpSeqPath, "seqB"), 'r') as fp:
            len_B=float("{:.8f}".format(len(fp.read())))
        
        avgSz=int((len_A+len_B)/2)
        os.system(f'rm {tmpSeqPath}/seqA {tmpSeqPath}/seqB')
        os.system(f'zcat {pairToMerge[0]} {pairToMerge[1]} | gto_fasta_to_seq | tr -d -c "ACGT" > {tmpSeqPath}/output.seq')
        os.system(f'{sodPath}SoD -u 7 a=1/1 d=1/100 -bs 50 -n {avgSz} -o {tmpSeqPath}/out.seq {tmpSeqPath}/output.seq')
        os.system(f'gto_fasta_from_seq < {tmpSeqPath}/out.seq > {tmpSeqPath}/{filename}')
        os.system(f'gzip {tmpSeqPath}/{filename}')
        os.replace(join(tmpSeqPath,filename + ".gz"),join(mergedPath,domain, filename + ".gz"))
        os.system(f'rm {tmpSeqPath}/output.seq {tmpSeqPath}/out.seq')
        counter+=1

def writeCSVLine(line):
    fileName = synthName+".csv"
    f = open(join(root, fileName), 'a')
    f.write(",".join(line))
    f.write("\n")
    f.close()

def getSequences():
    writeCSVLine(["Domain","DNA","AA","GC","L_DNA","L_AA"])
    for domain in listdir(mergedPath):
        print(f"Computing {domain}...")
        tmpPath = join(tmpDir,domain)
        if not path.exists(tmpPath):
            makedirs(tmpPath)

        for fileName in listdir(join(mergedPath, domain)):
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
            filePath = join(mergedPath, domain, fileName)
            if isfile(filePath):
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
                        print(f"WARNING: DNA of {name} was not computed! Maybe lack of free space in RAM!")
                        csvEntry["DNA"] = ""
                        
                os.system(f"AC -v -l 3 {tmpPath}/{name}.fna | sed '3,4d' | sed 'N;s/\\n/ /' > {tmpPath}/RESULTS_AC")
                with open(join(tmpPath, "RESULTS_AC"), 'r') as fp:
                    try:
                        csvEntry["AA"] = str(float(fp.read().split(" ")[15]))
                    except:
                        print(f"WARNING: AA of {name} was not computed! Maybe lack of free space in RAM!")
                        csvEntry["AA"] = ""
                os.system(f"rm {tmpPath}/*")
                writeCSVLine(csvEntry.values())



def _initialize():
    if not path.exists(root):
        makedirs(root)
    if not path.exists(tmpDir):
        makedirs(tmpDir)
    if not path.exists(tmpSeqPath):
        makedirs(tmpSeqPath)
    if not path.exists(mergedPath):
        makedirs(mergedPath)

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of SyntheticSequences/! There are relative paths defined in this code that need to be respected!")
    else:
        main()