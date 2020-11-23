#!/bin/bash

echo "Downloading in Data/Expression"
mkdir -p Data/Expression
cd Data/Expression
url="http://espinal3.inmegen.gob.mx:5270/genecorr/Expression"
wget -o down.log --cut-dirs=2 -r -nH -A '*-All.txt' $url
[[ $? -ne 0 ]] && echo "Download error: check down.log" && exit 15
echo "Succesfully download"
