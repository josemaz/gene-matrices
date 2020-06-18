#!/bin/bash

dirclean="Data/clean"

for t in Healthy Basal Her2 LumA LumB
do 
	julia julia/02-clusters.jl -c ${dirclean}/${t}-clean.tsv
done