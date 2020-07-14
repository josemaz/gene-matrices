#!/bin/bash

dirclean="Data/clean"

for t in Healthy Basal Her2 LumA LumB
do 
	echo "Using ${t} subtype"
	julia julia/02-clusters.jl -c ${dirclean}/${t}-clean.tsv
	julia julia/02-clusters-chrom.jl ${dirclean}/${t}-clean.tsv
done