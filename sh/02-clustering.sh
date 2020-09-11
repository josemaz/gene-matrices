#!/bin/bash

srcpath=$(dirname $(realpath $0))
. $srcpath/bash_colors.sh


dirclean="Data/clean"

t1=$SECONDS
for t in Healthy Basal Her2 LumA LumB
do 
	clr_bold clr_green "Using ${t} subtype"
	
	clr_bold clr_white "Cluster ALL genes"
	t2=$SECONDS
	julia julia/02-clusters.jl -c ${dirclean}/${t}-clean.tsv
	clr_bold clr_blue "[Cluster ALL genes] $(( ($SECONDS-t2) / 60)) min. elapsed"
	
	clr_bold clr_white "Cluster by chrom"
	t3=$SECONDS
	julia julia/02-clusters-chrom.jl ${dirclean}/${t}-clean.tsv
	clr_bold clr_blue "[Cluster by chrom] $(( ($SECONDS-t3) / 60)) min. elapsed"
done
clr_bold clr_blue "[Total time] $(( ($t1-a) / 60)) min. elapsed"
clr_bold clr_green "Successfully finished"