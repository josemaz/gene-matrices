# Image quality [https://www.shutterstock.com/blog/inches-to-pixels-resize-image-quality]

import pandas as pd
import numpy as np
import seaborn as sns
import sys, random, os
#from time import perf_counter
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from termcolor import colored, cprint
from pathlib import Path
from palettable.colorbrewer.qualitative import Dark2_7


np.set_printoptions(threshold=sys.maxsize)
logprint = lambda x: cprint(x, 'red', attrs=["bold"])

#return an adjancy matrix
def list2adj(dfl, start=1):
	M = np.empty([dfl.shape[0], dfl.shape[0]],dtype=int)
	cl = dfl.iloc[:,start].values.astype(int) # clusterid list set in column one.
	for i in range(len(cl)):
		clidsrc = cl[i]
		M[i,i] = clidsrc
		for j in range(i+1,len(cl)):	
			cliddst = cl[j]
			if clidsrc == cliddst:
				M[i,j] = clidsrc			
			else:
				M[i,j] = 0
			M[j,i] = M[i,j]
		if i % 1e3 == 0 and i != 0:
			print("list2adj iteration: ",i)
	return M

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""

    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:

    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    color_list[0] = [52/255, 58/255, 235/255, 1] # Pone el fondo azul
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def plotheat(data, oname, nc, label):

	mask_ut=np.triu(np.ones(A.shape)).astype(np.bool)
	# ax = sns.heatmap(A,xticklabels=False, yticklabels=False, cmap="Accent_r")
	cmaps = ['flag', 'prism', 'ocean', 'gist_earth', 'terrain', 'gist_stern',
            'gnuplot', 'gnuplot2', 'CMRmap', 'cubehelix', 'brg',
            'gist_rainbow', 'rainbow', 'jet', 'nipy_spectral', 'gist_ncar',
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn']
	cmaps = ['rainbow']
	for i in cmaps:
		name = oname + i + ".png"
		print(name)
		ax = sns.heatmap(data, mask=mask_ut, xticklabels=False, yticklabels=False, 
			cmap=discrete_cmap( nc, i))
		plt.title(label)
		# plt.ylabel('Gene start position')
		# plt.xlabel('Gene start position')
		plt.tight_layout()
		plt.savefig(name,dpi=300)
		plt.clf()

	# plt.show()



############################################################
###  MAIN
############################################################

#! python py/03-corr-heatmaps.py Data/clustering/Healthy-clusters.tsv
# fname = sys.argv[1] # tsv filename: gen,clusterid,cost,genestart,chromid

dirpng = "Plots/"
# healthyDir = "Data/Clustered/Healthy/Healthy-chr1-clusters.tsv"
healthyDir = "Data/Clustered/Healthy/"

chrs = np.arange(1, 22).tolist()
chrs.append('X') 

for chr in chrs:
	fname = healthyDir + "Healthy-chr" + str(chr) + "-clusters.tsv"
	label = fname.split('/')[-1].split('-')[0]
	cell = fname.split('/')[-1].split('-')[-2]
	d = dirpng + cell
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

	A = list2adj(healthy)


	# ! 2. PLOTTING
	oname = d + "/" + label + "-" + cell + "-"
	print("Plot name: " + oname)
	lab = 'Clusters for ' + label + ' - chr ' + str(chr)
	plotheat(A, oname, nclusters, lab)


	# #! 3. Save Order
	# #df.gname.to_csv("Data-Healthy-cluster-order.txt",index=False)

	subtypes = ["Basal","Her2","LumA","LumB"]
	for subtype in subtypes:
		logprint("Subtype: %s" % subtype)
		fname = "Data/Clustered/" + subtype + "/" + subtype + "-" + cell + "-clusters.tsv"
		df = pd.read_csv(fname, sep = ",")
		df = df.set_index('gname')
		df = df.loc[healthy.gname]
		A = list2adj(df, start=0)
		oname = d + "/" + subtype + "-" + cell + "-" #+ "-clusters.png"
		# print(oname)
		lab = 'Clusters for ' + subtype + ' - chr' + str(chr)
		plotheat(A,oname,nclusters,lab)




