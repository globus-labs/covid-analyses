#!/bin/bash

path=$1
FILE=$2
# create a temp file of just the smile column
cd $path


awk -F',' '{print $3 "\t" $2 "\t" $1}' ${FILE} > smiles_$FILE.smi
obabel smiles_$FILE.smi -O can_$FILE -ocan -e
awk '{print $3","$2","$1}' can_$FILE > new_$FILE

sed -i 's/[ \t]*$//' new_$FILE

rm smiles_$FILE.smi
rm can_$FILE

cd ..
