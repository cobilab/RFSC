from os import listdir, path, makedirs
from os.path import isfile, join
import random
from pprint import pprint
import sys
import shutil
import wget
import requests
from bs4 import BeautifulSoup

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
summariesPath = "./summaries"
dstPathOfOriginalSequences = "./original_sequences/"
locationOfDatabases = "../../References/"
multiColumnPos = 19

def main():
    _ensuringReproducibility()
    _initialize()
    selectedSequences = getSequences()
    downloadSequences(selectedSequences["from_url"], dstPathOfOriginalSequences)
    copySequences(selectedSequences["from_db"], locationOfDatabases, dstPathOfOriginalSequences)

def _ensuringReproducibility():
    random.seed(0)

def _initialize():
    if not path.exists(dstPathOfOriginalSequences):
        makedirs(dstPathOfOriginalSequences)

def getSequences():
    selectedSequences = {
        "from_url":{},
        "from_db":{}
    }
    onlyfiles = [f for f in listdir(summariesPath) if isfile(join(summariesPath, f)) and join(summariesPath, f).endswith("summary.txt")]
    for fileName in onlyfiles:
        with open(join(summariesPath, fileName), 'r') as fp:
            domain = fileName.split("_")[0]
            content = fp.readlines()
            lenOfFile = len(content)
            listOfValues = sorted(random.sample(range(0, lenOfFile), numberOfEntries[domain]))
            if "mitochondrion" in domain or "plastid" in domain:
                selectedSequences["from_db"][domain] = _getRandomEntries(listOfValues, content, True)
            else:
                selectedSequences["from_url"][domain] = _getRandomEntries(listOfValues, content)
    return selectedSequences

def _getRandomEntries(listOfValues, content, singleColumn=False):
    readsOfInterest = []
    for n in listOfValues:
        if singleColumn:
            readsOfInterest.append(content[n].replace("\n", ""))
        else:
            readsOfInterest.append(content[n].split("\t")[multiColumnPos])
    return readsOfInterest


def downloadSequences(selectedSequences, dst):
    for domain in selectedSequences:
        for entry in selectedSequences[domain]:
            r = requests.get(entry)
            soup = BeautifulSoup(r.text, 'html.parser')
            for link in soup.find_all('a'):
                if "genomic.fna.gz" in link.get('href') and "from_genomic" not in link.get('href'):
                    fLink = join(entry, link.get('href'))
                    print(f"Downloading {fLink}")
                    wget.download(fLink, out=dst)

def copySequences(selectedSequences, source, dst):
    for domain in selectedSequences:
        if domain == "plastid":
            mainPath = join(source, "NCBI-Plastid/NM-plastid/")
            _move(selectedSequences[domain], mainPath, dst)
        elif domain == "mitochondrion":
            mainPath = join(source, "NCBI-Mitochondrial/NM-mitochondrion/")
            _move(selectedSequences[domain], mainPath, dst)
        else:
            print(f"WARNING: Something is wrong! Domain {domain} not found")

def _move(selectedSequences, mainPath, dst):
    for entry in selectedSequences:
        path = join(mainPath, entry)
        fileDst = join(dst, entry)
        shutil.copyfile(path, fileDst)

if __name__ == "__main__":
    if "/" in sys.argv[0]:
        print("ERROR: Please run this script inside of src/SyntheticGenerator! There are relative paths defined in this code that need to be respected!")
    else:
        main()