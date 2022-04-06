from os import listdir
from os.path import isfile, join
import random
from pprint import pprint

#To ensure that the results are replicable
random.seed(0)

numberOfEntries = {
    "viral" : 20,
    "bacteria" : 20,
    "archaea" : 20,
    "fungi" : 20,
    "plant" : 10,
    "protozoa" : 10,
    "mitochondrion" : 20,
    "plastid" : 20
}
path = "./summaries"
multiColumnPos = 19

def main():
    selectedSequences = getSequences()
    pprint(selectedSequences)
    
def getSequences():
    selectedSequences = {}
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f)) and join(path, f).endswith("summary.txt")]
    for fileName in onlyfiles:
        with open(join(path, fileName), 'r') as fp:
            domain = fileName.split("_")[0]
            content = fp.readlines()
            lenOfFile = len(content)
            listOfValues = sorted(random.sample(range(0, lenOfFile), numberOfEntries[domain]))

            if "mitochondrion" in domain or "plastid" in domain:
                selectedSequences[domain] = _getRandomEntries(listOfValues, content, True)
            else:
                selectedSequences[domain] = _getRandomEntries(listOfValues, content)
    return selectedSequences

def _getRandomEntries(listOfValues, content, singleColumn=False):
    readsOfInterest = []
    for n in listOfValues:
        if singleColumn:
            readsOfInterest.append(content[n].replace("\n", ""))
        else:
            readsOfInterest.append(content[n].split("\t")[multiColumnPos])
    return readsOfInterest

if __name__ == "__main__":
    main()