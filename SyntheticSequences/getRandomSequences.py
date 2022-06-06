from os import listdir, path, makedirs
from os.path import isfile, join
import random
from pprint import pprint
import sys
import shutil
import wget
import requests
from bs4 import BeautifulSoup

#ensuring reproducibility
random.seed(0)

numberOfEntries = {
    "viral" : 7,
    "bacteria" : 7,
    "archaea" : 7,
    "fungi" : 7,
    "plant" : 7,
    "protozoa" : 7,
    "mitochondrion" : 7,
    "plastid" : 7
}
summariesPath = "../Mutations/summaries"
dstPathOfOriginalSequences = "./original_sequences/"
locationOfDatabases = "../References/"
multiColumnPos = 19

def main(split=False):
    selectedSequences = getSequences()
    if split:
        splitSequences(selectedSequences, dstPathOfOriginalSequences)
    else:
        _initialize()
        downloadSequences(selectedSequences["from_url"], dstPathOfOriginalSequences)
        copySequences(selectedSequences["from_db"], locationOfDatabases, dstPathOfOriginalSequences)

def _initialize():
    if (not path.exists(join(locationOfDatabases, "NCBI-Plastid/NM-plastid/"))) or len(listdir(join(locationOfDatabases, "NCBI-Plastid/NM-plastid/" ))) == 0:
        print("The plastid folder is empty, please download the database using the script:")
        print("./RFSC.sh  --download-ref-plastid ")
    if (not path.exists(join(locationOfDatabases, "NCBI-Mitochondrial/NM-mitochondrion/"))) or len(listdir(join(locationOfDatabases, "NCBI-Mitochondrial/NM-mitochondrion/" ))) == 0:
        print("The mitochondrion folder is empty, please download the database using the script:")
        print("./RFSC.sh  --download-ref-mitochondrial")
        sys.exit()
    if not path.exists(dstPathOfOriginalSequences):
        makedirs(dstPathOfOriginalSequences)
    
    

def splitSequences(selectedSequences, dstPathOfOriginalSequences):
    allFiles = [f for f in listdir(dstPathOfOriginalSequences) if isfile(join(dstPathOfOriginalSequences, f))]
    for x in selectedSequences:
        for domain in selectedSequences[x]:
            dst = join(dstPathOfOriginalSequences, domain)
            if not path.exists(dst):
                makedirs(dst)
            for entry in selectedSequences[x][domain]:
                fileName = entry
                if "/" in entry:
                    fileName = entry.split("/")[-1]
                fileToMove = list(filter(lambda x: fileName in x, allFiles))[0]
                shutil.move(join(dstPathOfOriginalSequences, fileToMove), dst)
    print("Done!")

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
        print("ERROR: Please run this script inside of Mutations/! There are relative paths defined in this code that need to be respected!")
    else:
        split = False
        if len(sys.argv) > 1:
            split = "-s" in sys.argv[1]
        main(split)