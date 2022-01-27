#!/bin/bash

#bash RFSC.sh --install

python3 -m pip install Conda/deps/wheel-0.22.0-py2.py3-none-any.whl
python3 -m pip install Conda/deps/keras-2.7.0-py2.py3-none-any.whl

wget https://files.pythonhosted.org/packages/52/4e/1d4186fc3cb6de68fe2572c7e148fabe70572608a46c7d2441ff74b56026/pandas-1.3.4-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl -O Conda/deps/pandas-1.3.4-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
python3 -m pip install Conda/deps/pandas-1.3.4-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl
rm Conda/deps/pandas-1.3.4-cp38-cp38-manylinux_2_17_x86_64.manylinux2014_x86_64.whl

wget https://files.pythonhosted.org/packages/36/66/eb5884679375385e316cee793db61c557e9b1f3b99b8559c88ac14d7e561/numpy-1.21.4-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl -O Conda/deps/numpy-1.21.4-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl
python3 -m pip install Conda/deps/numpy-1.21.4-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl
rm Conda/deps/numpy-1.21.4-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl

wget https://files.pythonhosted.org/packages/ae/bf/e1ce8638e6c5234632c4d87635ba23f4b4b837ba8c85ba8f8262569cdeb1/scikit_learn-1.0.1-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl -O Conda/deps/scikit_learn-1.0.1-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl
python3 -m pip install Conda/deps/scikit_learn-1.0.1-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl
rm Conda/deps/scikit_learn-1.0.1-cp38-cp38-manylinux_2_12_x86_64.manylinux2010_x86_64.whl

wget https://files.pythonhosted.org/packages/72/8a/033b584f8dd863c07aa8877c2dd231777de0bb0b1338f4ac6a81999980ee/tensorflow-2.7.0-cp38-cp38-manylinux2010_x86_64.whl -O Conda/deps/tensorflow-2.7.0-cp38-cp38-manylinux2010_x86_64.whl
python -m pip install Conda/deps/tensorflow-2.7.0-cp38-cp38-manylinux2010_x86_64.whl
rm Conda/deps/tensorflow-2.7.0-cp38-cp38-manylinux2010_x86_64.whl

wget https://files.pythonhosted.org/packages/a7/c9/4968ca0434c313aed71fc4cc2339aa8844482d5eefdcc8989c985a19ea2e/xgboost-1.5.0-py3-none-manylinux2014_x86_64.whl -O Conda/deps/xgboost-1.5.0-py3-none-manylinux2014_x86_64.whl
python3 -m pip install Conda/deps/xgboost-1.5.0-py3-none-manylinux2014_x86_64.whl
rm Conda/deps/xgboost-1.5.0-py3-none-manylinux2014_x86_64.whl

wget http://cab.spbu.ru/files/release3.15.3/SPAdes-3.15.3-Linux.tar.gz
tar -xzf SPAdes-3.15.3-Linux.tar.gz
cp SPAdes-3.15.3-Linux/bin/cds-mapping-stats $PREFIX/bin/cds-mapping-stats
cp SPAdes-3.15.3-Linux/bin/cds-subgraphs $PREFIX/bin/cds-subgraphs
cp SPAdes-3.15.3-Linux/bin/coronaspades.py $PREFIX/bin/coronaspades.py
cp SPAdes-3.15.3-Linux/bin/mag-improve $PREFIX/bin/mag-improve
cp SPAdes-3.15.3-Linux/bin/metaplasmidspades.py $PREFIX/bin/metaplasmidspades.py
cp SPAdes-3.15.3-Linux/bin/metaspades.py $PREFIX/bin/metaspades.py
cp SPAdes-3.15.3-Linux/bin/metaviralspades.py $PREFIX/bin/metaviralspades.py
cp SPAdes-3.15.3-Linux/bin/plasmidspades.py $PREFIX/bin/plasmidspades.py
cp SPAdes-3.15.3-Linux/bin/rnaspades.py $PREFIX/bin/rnaspades.py
cp SPAdes-3.15.3-Linux/bin/rnaviralspades.py $PREFIX/bin/rnaviralspades.py
cp SPAdes-3.15.3-Linux/bin/spades-bwa $PREFIX/bin/spades-bwa
cp SPAdes-3.15.3-Linux/bin/spades-convert-bin-to-fasta $PREFIX/bin/spades-convert-bin-to-fasta
cp SPAdes-3.15.3-Linux/bin/spades-core $PREFIX/bin/spades-core
cp SPAdes-3.15.3-Linux/bin/spades-corrector-core $PREFIX/bin/spades-corrector-core
cp SPAdes-3.15.3-Linux/bin/spades-gbuilder $PREFIX/bin/spades-gbuilder
cp SPAdes-3.15.3-Linux/bin/spades-gmapper $PREFIX/bin/spades-gmapper
cp SPAdes-3.15.3-Linux/bin/spades-gsimplifier $PREFIX/bin/spades-gsimplifier
cp SPAdes-3.15.3-Linux/bin/spades-hammer $PREFIX/bin/spades-hammer
cp SPAdes-3.15.3-Linux/bin/spades_init.py $PREFIX/bin/spades_init.py
cp SPAdes-3.15.3-Linux/bin/spades-ionhammer $PREFIX/bin/spades-ionhammer
cp SPAdes-3.15.3-Linux/bin/spades-kmercount $PREFIX/bin/spades-kmercount
cp SPAdes-3.15.3-Linux/bin/spades-kmer-estimating $PREFIX/bin/spades-kmer-estimating
cp SPAdes-3.15.3-Linux/bin/spades.py $PREFIX/bin/spades.py
cp SPAdes-3.15.3-Linux/bin/spades-read-filter $PREFIX/bin/spades-read-filter
cp SPAdes-3.15.3-Linux/bin/spades-truseq-scfcorrection $PREFIX/bin/spades-truseq-scfcorrection
cp SPAdes-3.15.3-Linux/bin/spaligner $PREFIX/bin/spaligner
cp SPAdes-3.15.3-Linux/bin/truspades.py $PREFIX/bin/truspades.py
rm -r SPAdes-3.15.3-Linux

cp Conda/deps/blastn $PREFIX/bin/blastn

bash RFSC.sh --conda

cp RFSC.sh $PREFIX/bin/
