
# gene-matrices
This repository contains code and supplementary materials for paper named "xxx".  Alfredo Gonzalez, Jose Maria Zamora-Fuentes, Enrique Hernandez-Lemus, Jesus Espinal-Enriquez


## Install enviroment

Install miniconda:

`$ bash install-miniconda.sh`

Install julia:

`$ bash install-julia.sh`



## Bioinformatics process


### 0.5 - Download data

`$ bash download-data.sh`


### 01 - Clean subtype data

`$ bash 01-cleandata.sh`

Outputs are saved on:

`$ ls Data/clean/  `


### 02 - Clustering 

This step can be slow.

`$ bash 02-clustering.sh`

Execute in backgorund to best performance:

`$ bash 02-clustering.sh &> salida.log &`

Clustering Outputs are saved on:

`$ ls Data/Clustered  `

Clustering by independent chromosome was named with *chr* label, and clustering by all chromosomes was named as *all* label in each subtype directory.

Pearson triangular matrix in list format are saved on:

`$ ls Data/Pearson/  `


### 03 - Pearson correlation

`$ python py/03-plot-pearson.py Data/clean/Healthy-clean.tsv`

Plots are saved on:

`$ ls Plots/Healthy/    # for example`