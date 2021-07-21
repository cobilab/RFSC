#!/bin/bash

# Script used to output the predictions percentages presented
# in a file previously generated (i.e. ./Tests/deepTest_(...).sh)
#
# Outputs the solution in the terminal with all the percentages
# of each domain
#
# Run:      ./Tests/predictions.sh FILE
# Example:  ./Tests/predictions.sh Prediction_3_7.txt

INPUT_FILE=$1;

VIRUS=0;
BACTERIA=0;
ARCHAEA=0;
FUNGI=0;
PLANT=0;
PROTOZOA=0;
MITO=0;
PLASTID=0;

TOTAL=0;

file_name=$(basename $INPUT_FILE);

VIRUS=`grep -o 'Virus\|Viral' $INPUT_FILE | wc -l`
BACTERIA=`grep -o 'Bacteria' $INPUT_FILE | wc -l`
ARCHAEA=`grep -o 'Archaea' $INPUT_FILE | wc -l`
FUNGI=`grep -o 'Fungi' $INPUT_FILE | wc -l`
PLANT=`grep -o 'Plant' $INPUT_FILE | wc -l`
PROTOZOA=`grep -o 'Protozoa' $INPUT_FILE | wc -l`
MITO=`grep -o 'Mitochondrial' $INPUT_FILE | wc -l`
PLASTID=`grep -o 'Plastid' $INPUT_FILE | wc -l`

TOTAL=$(( $VIRUS + $BACTERIA + $ARCHAEA + $FUNGI + $PLANT + $PROTOZOA + $MITO + $PLASTID ))


VIRUS=`echo "scale=3 ; $VIRUS / $TOTAL * 100" | bc`
BACTERIA=`echo "scale=3 ; $BACTERIA / $TOTAL * 100" | bc`
ARCHAEA=`echo "scale=3 ; $ARCHAEA / $TOTAL * 100" | bc`
FUNGI=`echo "scale=3 ; $FUNGI / $TOTAL * 100" | bc`
PLANT=`echo "scale=3 ; $PLANT / $TOTAL * 100" | bc`
PROTOZOA=`echo "scale=3 ; $PROTOZOA / $TOTAL * 100" | bc`
MITO=`echo "scale=3 ; $MITO / $TOTAL * 100" | bc`
PLASTID=`echo "scale=3 ; $PLASTID / $TOTAL * 100" | bc`

echo "                                                                    "
echo " Prediction results from $file_name: *****************              "
echo "                                                                    "
echo " *** Viral Percentage: .................................. (%) $VIRUS "
echo " *** Bacteria Percentage: ............................... (%) $BACTERIA "
echo " *** Archaea Percentage: ................................ (%) $ARCHAEA "
echo " *** Fungi Percentage: .................................. (%) $FUNGI "
echo " *** Plant Percentage: .................................. (%) $PLANT "
echo " *** Protozoa Percentage: ............................... (%) $PROTOZOA "
echo " *** Mitochondrial Percentage: .......................... (%) $MITO "
echo " *** Plastid Percentage: ................................ (%) $PLASTID "
echo "                                                                    "
echo " *****************************************************************  "
echo "                                                                    "