# gene-matrices
This repository contains code and supplementary materials for paper named "xxx".  Alfredo Gonzalez, Jose Maria Zamora-Fuentes, Enrique Hernandez-Lemus, Jesus Espinal-Enriquez

## Install enviroment

Install miniconda:

`$ bash install-miniconda.sh`

Install julia:

`$ bash install-julia.sh`

Install data directory:

`$ bash install-data.sh`

## 01 - Clean subtype data

`$ bash 01-clean.sh`

## 02 - Clustering and Pearson correlation

`$ julia julia julia/02-clusters.jl -c Data/clean/Healthy-clean.tsv`

`$ bash 02-CorrAndPearson.sh`
