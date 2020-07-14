#!/bin/bash

dirclean="Data/clean"

for t in Healthy Basal Her2 LumA LumB
do 
	echo "Using ${t} subtype"
	python py/03-plot-pearson.py Data/clean/${t}-clean.tsv
done