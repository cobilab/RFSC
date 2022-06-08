<p align="center"> <img src="Logo.png"> </p>


[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)


RFSC is a Reference-Free Sequence Classification Tool that using machine learning classifiers relies on an ensemble of experts in order to provide efficient classification in metagenomic contexts.


## <b>Installation</b>

### <b>Using Docker</b>

```sh
git clone https://github.com/cobilab/RFSC
cd RFSC
chmod +x RFSC.sh 
chmod +777 src/*.sh
docker-compose build
docker-compose up -d && docker exec -it rfsc bash && docker-compose down
./RFSC.sh --install #install tools
```

### <b>Download NCBI Reference Databases</b>

```sh
./RFSC.sh --build-ref-virus --build-ref-bacteria --build-ref-archaea --build-ref-protozoa \ --build-ref-fungi --build-ref-plant --build-ref-mitochondrial --build-ref-plastid
```

or 

```sh
./RFSC.sh -dviral -dbact -darch -dprot -dfung -dplan -dmito -dplas
```

# Global Results

### <b>Real Sequence Classification </b>
Obtain classification report of KNN, GNB and XGBoost. 

```sh
./RFSC.sh -runAll #classification report table
./RFSC.sh -runAll F1Score # Weighted-averaged F1-score
./RFSC.sh -runAll Accuracy # Average Accuracy
```

### <b>Generate mutated data, and perform classification</b>

#### To gathers a small set of sequences from the 8 domains, run the script:
```sh
./RFSC.sh -mget  
```

#### To compute features for mutated sequences, run the script:
```sh
./RFSC.sh -cfem
```

#### To perform classification for all mutated sequences, run the script:
```sh
./RFSC.sh -cclm 
```

### <b>Synthetic Sequence Generation and Classification</b>

##### To gathers a small set of sequences from the 8 domains, run the script:
```sh
./RFSC.sh -sget
```

#### To create the synthetic hybrid sequences and compute their features, run the script:

```sh
./RFSC.sh -cfes
```

#### To perform classification of the synthetic hybrid sequences and obtain classification report of KNN, GNB and XGBoost, run the script:

```sh
./RFSC.sh -ccls
```

#### To perform classification of the synthetic sequences using Kraken2, run the script:
#####  (only for comparison purposes, requires Kraken2 installation)
##### You should download the Kraken2 database at: https://benlangmead.github.io/aws-indexes/k2 
##### To obtain the same results, use the Standard database containing "archaea, bacteria, viral, plasmid, human1, UniVec_Core" created at 5/17/2021, with 38.6GB.
```sh
./RFSC.sh -ckra
```

# Running Examples

##### ✨ Generate a synthetic sequence and subsequently proceed to a Reference-Free Reconstruction of the same:
&nbsp;
```sh
./RFSC.sh --clean y
./RFSC.sh --threads 8 --gen-adapters
./RFSC.sh --efetch-fasta 155971 Input_Data/EntrezGenomes 
./RFSC.sh --efetch-fasta EF491856.1 Input_Data/EntrezGenomes 
./RFSC.sh --efetch-fasta MT682520 Input_Data/EntrezGenomes
./RFSC.sh -synt Input_Data/EntrezGenomes/155971.fna Input_Data/EntrezGenomes/EF491856.1.fna Input_Data/EntrezGenomes/MT682520.fna
./RFSC.sh -trim TT PE --run-de-novo
```

##### ✨ Reference-Based Classification, usign FALCON-meta:
###### (If the reference databases have already been built and the Reference Free Reconstruction stage is finished)
&nbsp;
```sh
./RFSC.sh --threads 8 --set-len-cov 100 3 --set-threshold-max-min 70 1 --run-falcon SO Viral
```


##### ✨ Reference-Free Classification, using XBoost 
&nbsp;
```sh
./RFSC.sh --threads 8 --efetch-fasta 155971 RefFree
./RFSC.sh --run-xgboost
```

# System Requirements

Laptop computer running Linux Ubuntu (for example, 18.04 LTS or higher) with GCC (https://gcc.gnu.org), Conda (https://docs.conda.io) and CMake (https://cmake.org) installed. The hardware must contain at least 32 GB of RAM, and a 800 GB disk. In the case of the this, if the database is not re-built, it is only needed near 10 GB of space.

## Tools Integrated in RFSC

| Tool | URL |
| ------ | ------ |
| AC | [https://github.com/cobilab/ac][PlGa] |
| Blastn | [https://blast.ncbi.nlm.nih.gov/Blast.cgi][PlOd] |
| Cryfa | [https://github.com/cobilab/cryfa][PlGd] |
| Entrez | [https://www.ncbi.nlm.nih.gov/genome][PlMe] |
| FALCON-meta | [https://github.com/cobilab/falcon][PlGa] |
| FASTP | [https://github.com/OpenGene/fastp][PlGh] |
| GeCo3 | [https://github.com/cobilab/geco3][PlMe] |
| GTO | [https://cobilab.github.io/gto/][PlOd] |
| metaSPAdes | [https://cab.spbu.ru/software/meta-spades/][PlGd] |
| ORFfinder | [https://www.ncbi.nlm.nih.gov/orffinder/][PlMe] |
| ORFM | [https://github.com/wwood/OrfM][PlGa] |
| SoD | [https://github.com/pratas/SoD.git][PlDb] |
| Trimmomatic | [http://www.usadellab.org/cms/?page=trimmomatic][PlDb] |
## License

GNU GPL

**✨Developed to make a change!✨**

[//]: # (These are reference links used in the body of this note and get stripped out when the markdown processor does its job. There is no need to format nicely because it shouldn't be seen. Thanks SO - http://stackoverflow.com/questions/4823468/store-comments-in-markdown-syntax)

   [dill]: <https://github.com/joemccann/dillinger>
   [git-repo-url]: <https://github.com/joemccann/dillinger.git>
   [john gruber]: <http://daringfireball.net>
   [df1]: <http://daringfireball.net/projects/markdown/>
   [markdown-it]: <https://github.com/markdown-it/markdown-it>
   [Ace Editor]: <http://ace.ajax.org>
   [node.js]: <http://nodejs.org>
   [Twitter Bootstrap]: <http://twitter.github.com/bootstrap/>
   [jQuery]: <http://jquery.com>
   [@tjholowaychuk]: <http://twitter.com/tjholowaychuk>
   [express]: <http://expressjs.com>
   [AngularJS]: <http://angularjs.org>
   [Gulp]: <http://gulpjs.com>

   [PlDb]: <https://github.com/joemccann/dillinger/tree/master/plugins/dropbox/README.md>
   [PlGh]: <https://github.com/joemccann/dillinger/tree/master/plugins/github/README.md>
   [PlGd]: <https://github.com/joemccann/dillinger/tree/master/plugins/googledrive/README.md>
   [PlOd]: <https://github.com/joemccann/dillinger/tree/master/plugins/onedrive/README.md>
   [PlMe]: <https://github.com/joemccann/dillinger/tree/master/plugins/medium/README.md>
   [PlGa]: <https://github.com/RahulHP/dillinger/blob/master/plugins/googleanalytics/README.md>
