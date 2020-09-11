#!/bin/bash

condadir="${HOME}/bioconda"
cd /tmp
pkg="Miniconda3-latest-Linux-x86_64.sh"
pkg="Miniconda3-py37_4.8.3-Linux-x86_64.sh"
wget "https://repo.anaconda.com/miniconda/${pkg}"
bash ${pkg} -p $condadir -b -s
rm -rf /tmp/${pkg}
. $condadir/bin/activate
conda install -y numpy pandas seaborn palettable termcolor
easy_install trash-cli
pip install python-igraph
pip install gprofiler-official
echo ". $condadir/bin/activate" >> ~/.bashrc
