#!/bin/bash
#
# ==============================================================================
#
#   IT DOWNLOADS AND COMPRESSES THE COMPLETE GENOMIC 
#   SEQUENCES FROM THE NCBI DATABASE.
#
#   THE OUTPUT FILES ARE:
#
#   DB-viral.fa.gz                 (VIRAL)
#   DB-bacteria.fa.gz              (BACTERIA)
#   DB-archaea.fa.gz               (ARCHAEA)
#   DB-protozoa.fa.gz              (PROTOZOA)
#   DB-fungi.fa.gz                 (FUNGI)
#   DB-plant.fa.gz                 (PLANT)
#   DB-investebrate.fa.gz          (INVERTEBRATE)
#   DB-vertebrate_mammalian.fa.gz  (VERTEBRATE MAMMALIAN)
#   DB-vertebrate_other.fa.gz      (VERTEBRATE OTER)
#   DB-mitochondrion.fa.gz         (MITOCHONDRIAL)
#   DB-plastid.fa.gz               (PLASTID)
#
# ==============================================================================
#
THREADS=4;
VIRAL=0;
BACTERIA=0;
ARCHAEA=0;
PROTOZOA=0;
FUNGI=0;
PLANT=0;
INVERTEBRATE=0;
VERTEBRATE_MAMMALIAN=0;
VERTEBRATE_OTHER=0;
MITOCHONDRION=0;
PLASTID=0;
HELP=0;
#
# ==============================================================================
#
if [ "$#" -eq 0 ];
  then
  HELP=1;
  fi
#
POSITIONAL=();
#
while [[ $# -gt 0 ]]
  do
  i="$1";
  case $i in
    -h|--help)                                          HELP=1; shift  ;;
    -t|--threads)               THREADS="$2";           HELP=0; shift 2;;
    -vi|--viral)                VIRAL=1;                HELP=0; shift  ;;
    -ba|--bacteria)             BACTERIA=1;             HELP=0; shift  ;;
    -ar|--archaea)              ARCHAEA=1;              HELP=0; shift  ;;
    -pr|--protozoa)             PROTOZOA=1;             HELP=0; shift  ;;
    -fu|--fungi)                FUNGI=1;                HELP=0; shift  ;;
    -pl|--plant)                PLANT=1;                HELP=0; shift  ;;
    -in|--invertebrate)         INVERTEBRATE=1;         HELP=0; shift  ;;
    -vm|--vertebrate_mammalian) VERTEBRATE_MAMMALIAN=1; HELP=0; shift  ;;
    -vo|--vertebrate_other)     VERTEBRATE_OTHER=1;     HELP=0; shift  ;;
    -mi|--mitochondrion )       MITOCHONDRION=1;        HELP=0; shift  ;;
    -ps|--plastid)              PLASTID=1;              HELP=0; shift  ;;
    *) echo "Invalid arg "$1 ; exit 1;
    ;;
  esac
done
set -- "${POSITIONAL[@]}" # restore positional parameters
#
# ==============================================================================
# HELP
#
if [[ "$HELP" -eq "1" ]]
  then
  echo "                                                                    "
  echo "Usage: BUILD_DB.sh [options]                              "
  echo "                                                                    "
  echo "It downloads and compresses the choosen NCBI database:              "
  echo "                                                                    "
  echo "   -h,   --help                   Show this help message,           "
  echo "   -t,   --threads <INT>          Number of threads (DEF: 4),       "
  echo "                                                                    "
  echo "   -vi,  --viral                  Builds viral DB,                  "
  echo "   -ba,  --bacteria               Builds bacteria DB,               "
  echo "   -ar,  --archaea                Builds archaea DB,                "
  echo "   -pr,  --protozoa               Builds protozoa DB,               "
  echo "   -fu,  --fungi                  Builds fungi DB,                  "
  echo "   -pl,  --plant                  Builds plant DB,                  "
  echo "   -in,  --invertebrate           Builds invertebrate DB,           "
  echo "   -vm,  --vertebrate_mammalian   Builds vertebrate mammalian DB,   "
  echo "   -vo,  --vertebrate_other       Builds vertebrate other DB,       "
  echo "   -mi,  --mitochondrion          Builds mitochondrion DB,          "
  echo "   -ps,  --plastid                Builds plastid DB.                "
  echo "                                                                    "
  echo "Example: BUILD_DB.sh --viral --fungi --plant              "
  echo "                                                                    "
  exit 1
  fi
#
# ==============================================================================
#
function BUILD_DB () {
#
echo -e "\033[1;34m[RFSC] \033[1;32m Downloading $1 DB using $2 threads ... \033[0m";
#
rm -f $1-url_download.txt;
mkdir -p NM-$1;
#
curl 'ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/'$1'/assembly_summary.txt' | \
awk '{FS="\t"} !/^#/ {print $20} ' | \
sed -r 's|(ftp://ftp.ncbi.nlm.nih.gov/genomes/all/.+/)(GCF_.+)|\1\2/\2_genomic.fna.gz|' > $1-url_download.txt
### awk -F '\t' '{if($12=="Complete Genome"||$12=="Chromosome") print $20}' assembly_summary.txt
#
cat $1-url_download.txt | xargs -n 1 -P $2 wget --directory-prefix=NM-$1 -q
#
rm -f DB-$1.fa $1-url_download.txt;
#XXX: VERIFY ERRORS (> in the middle of the DNA and Gzip errors...)
for f in NM-$1/*.fna.gz; do zcat "$f" >> DB-$1.fa; done
#rm -fr NM-$1;
#
gzip -f DB-$1.fa;
#
echo -e "\033[1;34m[RFSC] \033[1;32m Building has been successful! \033[0m";
#
}
#
# ==============================================================================
#
if [[ "$VIRAL" -eq "1" ]];
  then
  BUILD_DB "viral" "$THREADS";
  fi
#
if [[ "$BACTERIA" -eq "1" ]];
  then
  BUILD_DB "bacteria" "$THREADS";
  fi
#
if [[ "$ARCHAEA" -eq "1" ]];
  then
  BUILD_DB "archaea" "$THREADS";
  fi
#
if [[ "$PROTOZOA" -eq "1" ]];
  then
  BUILD_DB "protozoa" "$THREADS";
  fi
#
if [[ "$FUNGI" -eq "1" ]];
  then
  BUILD_DB "fungi" "$THREADS";
  fi
#
if [[ "$PLANT" -eq "1" ]];
  then
  BUILD_DB "plant" "$THREADS";
  fi
#
if [[ "$INVERTEBRATE" -eq "1" ]];
  then
  BUILD_DB "invertebrate" "$THREADS";
  fi
#
if [[ "$VERTEBRATE_MAMMALIAN" -eq "1" ]];
  then
  BUILD_DB "vertebrate_mammalian" "$THREADS";
  fi
#
if [[ "$VERTEBRATE_OTHER" -eq "1" ]];
  then
  BUILD_DB "vertebrate_other" "$THREADS";
  fi
#
if [[ "$MITOCHONDRION" -eq "1" ]];
  then
  echo -e "\033[1;34m[RFSC] \033[1;32m Downloading mitochondrion DB ...\033[0m";
  MT_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion";
  rm -f MT1.fna.gz MT2.fna.gz DB-mitochondrion.fa.gz
  wget $MT_URL/mitochondrion.1.1.genomic.fna.gz -O MT1.fna.gz
  wget $MT_URL/mitochondrion.2.1.genomic.fna.gz -O MT2.fna.gz
  zcat MT1.fna.gz MT2.fna.gz | gzip -f -5 > DB-mitochondrion.fa.gz
  rm -f MT1.fna.gz MT2.fna.gz
  echo -e "\033[1;34m[RFSC] \033[1;32m Building has been successful! \033[0m";
  fi
#
if [[ "$PLASTID" -eq "1" ]];
  then
  echo -e "\033[1;34m[RFSC] \033[1;32m Downloading plastid DB ...\033[0m";
  PL_URL="ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plastid";
  rm -f PL1.fna.gz PL2.fna.gz DB-plastid.fa.gz
  wget $PL_URL/plastid.1.1.genomic.fna.gz -O PL1.fna.gz
  wget $PL_URL/plastid.2.1.genomic.fna.gz -O PL2.fna.gz
  wget $PL_URL/plastid.3.1.genomic.fna.gz -O PL3.fna.gz
  wget $PL_URL/plastid.4.1.genomic.fna.gz -O PL4.fna.gz
  zcat PL1.fna.gz PL2.fna.gz PL3.fna.gz PL4.fna.gz \
  | gzip -f -5 > DB-plastid.fa.gz
  rm -f PL1.fna.gz PL2.fna.gz PL3.fna.gz PL4.fna.gz
  echo -e "\033[1;34m[RFSC] \033[1;32m Building has been successful! \033[0m";
  fi
#
# ==============================================================================
#