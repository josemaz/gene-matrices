#!/bin/bash

condadir="${HOME}/bioconda"
cd /tmp
wget "https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh"
bash Miniconda3-latest-Linux-x86_64.sh -p $condadir -b -s
rm -rf /tmp/Miniconda3-latest-Linux-x86_64.sh
. $condadir/bin/activate
conda install -y numpy pandas seaborn
easy_install trash-cli
pip install python-igraph
pip install gprofiler-official
echo ". $condadir/bin/activate" >> ~/.bashrc
