import numpy as np
from termcolor import colored, cprint
from pathlib import Path
import pandas as pd
import sys

logprint = lambda x: cprint(x, 'red', attrs=["bold"])

dirPlots = "Plots/"
dirHealthy = "Data/Clustered/Healthy/"

chrs = np.arange(1, 22).tolist()
chrs.append('X') 

for chr in chrs:
	fname = dirHealthy + "Healthy-chr" + str(chr) + "-clusters.tsv"
	label = fname.split('/')[-1].split('-')[0]
	cell = fname.split('/')[-1].split('-')[-2]
	d = dirPlots + cell
	logprint("Plot dir: %s" % d)
	logprint("Plot label: %s" % label)
	logprint("Plot chr name: %s" % cell)
	fp = Path(d)
	fp.mkdir(parents=True, exist_ok=True)


	#! 1. PREPARING DATA
	healthy = pd.read_csv(fname, sep = ",")
	nclusters = healthy["clusterid"].max()

	# healthy.sort_values(by=['chromosome','gstart']) # Order by cromosome and then gene start
	healthy = healthy.sort_values(by=['chromosome','clusterid']) # Order by cromosome and clusterid
	print("Number of Clusters: %s" % nclusters)
	print("Number of Genes: %s" % healthy.shape[0]) # Print rows

	sys.exit(15)
