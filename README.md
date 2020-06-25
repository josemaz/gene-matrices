# gene-matrices
This repository contains code and supplementary materials for paper named "xxx".  Alfredo Gonzalez, Jose Maria Zamora-Fuentes, Enrique Hernandez-Lemus, Jesus Espinal-Enriquez


## Install enviroment

Install miniconda:

`$ bash install-miniconda.sh`

Install julia:

`$ bash install-julia.sh`

Install data directory:

`$ bash install-data.sh`


## Bioinformatics process

# 01 - Clean subtype data

`$ bash 01-clean.sh`

# 02 - Clustering 

`$ bash 02-clustering.sh`

Execute in backgorund to best performance:

`$ bash 02-clustering.sh &> salida.log &`

# Pearson correlation

`$ python py/03-plot-pearson.py Data/clean/Healthy-clean.tsv`

Plots are saved on:

`$ ls Plots/Healthy/    # for example`