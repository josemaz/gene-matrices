
# gene-matrices
This repository contains code and supplementary materials for paper named "xxx".  Alfredo Gonzalez, Jose Maria Zamora-Fuentes, Enrique Hernandez-Lemus, Jesus Espinal-Enriquez



## Pre-requisites: nstall enviroment

Install miniconda:

`$ bash sh/install-miniconda.sh`

Install julia:

`$ bash sh/install-julia.sh`



## Bioinformatic process


### 0.5 - Download data

We shared expression matrices as public data in:

http://espinal3.inmegen.gob.mx:5270/genecorr/subtipos-originales/

to save locally on `Data/Expression`:

`$ bash sh/download-data.sh`


### 01 - Clean subtype data

`$ bash sh/01-cleandata.sh`

Outputs are saved on:

`$ ls Data/clean/`


### 02 - Clustering 

This step can be slow.

`$ bash sh/02-clustering.sh`

Execute in backgorund to best performance:

`$ bash sh/02-clustering.sh &> salida.log &`

Clustering Outputs are saved on:

`$ ls Data/Clustered  `

Clustering by independent chromosome was named with *chr* label, and clustering by all chromosomes was named as *all* label in each subtype directory.

Pearson triangular matrices are in list format are saved on:

`$ ls Data/Pearson/  `

**Optional:** If you don't have computing hardware. You can download our results with

`$ bash sh/download-pearson.sh`

**Warnning:** Pearson correlations can use 30GB of space.



## Visualization

### 03 - Pearson correlation  

To get heatmaps of pearson corelation for every subtype: 

`$ bash sh/03-plot-corr.sh`

Plots are saved on:

`$ ls Plots/Healthy/    # for example`

### 04 - Plot clusters

`$ python py/04-plot-clusters.py`

Plots are saved by chromosome as:

`$ ls Plots/chr1/Healthy-chr1-rainbow.png    # for example`
`$ ls Plots/chr1/Basal-chr1-rainbow.png    # for example`

### 05 - Plot correlations and cumulative distribution

`$ python py/05-plot-gstart.py`

Plots are saved by chromosome as:

`$ ls Plots/chr1/Healthy-chr1-gstart-cum.png`
`$ ls Plots/chr1/Healthy-chr1-gstart-heat.png`

Files with `cum` suffix are cummilative plots of elements in every cluster. Moreover `heat` suffix are heatmap plots by cluster. The colors are the same for each cluster in `cum` and `heat`.






