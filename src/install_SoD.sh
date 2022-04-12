#!/bin/bash
# =============================================
# Instalation of SoD
# =============================================
rm -rf SoD
mkdir SoD
cd SoD

git clone https://github.com/pratas/SoD.git;
pwd
cd SoD/src/
make 
cd ../../
mv  SoD/src/SoD ./sod
rm -rf SoD
mv sod SoD
chmod +x SoD
