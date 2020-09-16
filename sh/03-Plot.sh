#!/bin/bash

srcpath=$(dirname $(realpath $0))
. $srcpath/bash_colors.sh

dirclean="Data/clean"

for t in Healthy Basal Her2 LumA LumB
do 
	clr_bold clr_green "Using ${t} subtype"
	python py/03-plot-pearson.py Data/clean/${t}-clean.tsv
done

# Plots of clusters and cumulative distributions
# Healthy colors vs subtypes
python py/04-plot-clusters.py