<p align="center"> <img src="Logo.png"> </p>


[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](LICENSE)


RFSC is a Reference-Free Sequence Classification Tool that using machine learning classifiers relies on an ensemble of experts in order to provide efficient classification in metagenomic contexts.


## Instalation

```sh
git clone https://github.com/cobilab/RFSC
cd RFSC
./RFSC.sh --install
```

## Using Docker

```sh
git clone https://github.com/cobilab/RFSC
cd RFSC
docker-compose build
docker-compose up -d && docker exec -it rfsc bash && docker-compose down
```

## Build NCBI Reference Databases

```sh
./RFSC.sh --build-ref-virus --build-ref-bacteria --build-ref-archaea --build-ref-protozoa \ --build-ref-fungi --build-ref-plant --build-ref-mitochondrial --build-ref-plastid
```

## Running Examples

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

## System Requirements

Laptop computer running Linux Ubuntu (for example, 18.04 LTS or higher) with GCC (https://gcc.gnu.org), Conda (https://docs.conda.io) and CMake (https://cmake.org) installed. The hardware must contain at least 8 GB of RAM, and a 800 GB disk. In the case of the this, if the database is not re-built, it is only needed near 10 GB of space.

## Tools Integrated in RFSC

| Tool | URL |
| ------ | ------ |
| Trimmomatic | [http://www.usadellab.org/cms/?page=trimmomatic][PlDb] |
| FASTP | [https://github.com/OpenGene/fastp][PlGh] |
| metaSPAdes | [https://cab.spbu.ru/software/meta-spades/][PlGd] |
| GTO | [https://cobilab.github.io/gto/][PlOd] |
| Entrez | [https://www.ncbi.nlm.nih.gov/genome][PlMe] |
| FALCON-meta | [https://github.com/cobilab/falcon][PlGa] |
| Cryfa | [https://github.com/cobilab/cryfa][PlGd] |
| Blastn | [https://blast.ncbi.nlm.nih.gov/Blast.cgi][PlOd] |
| ORFfinder | [https://www.ncbi.nlm.nih.gov/orffinder/][PlMe] |
| ORFM | [https://github.com/wwood/OrfM][PlGa] |
| GeCo3 | [https://github.com/cobilab/geco3][PlMe] |
| AC | [https://github.com/cobilab/ac][PlGa] |

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
