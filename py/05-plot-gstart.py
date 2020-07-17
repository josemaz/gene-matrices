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
    color_list[0] = [1, 1, 1, 1] # Pone el fondo azul
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)


def cmapfixed(ncolors):
	fixedcolors =10
	base = plt.cm.get_cmap('rainbow') #paired
	color_list = base(np.linspace(0, 1, fixedcolors))
	# color_list = np.insert(color_list,0,[1, 1, 1, 1], axis=0) # Pone el fondo blanco
	color_list = np.insert(color_list,0,[21/255, 25/255, 235/255, 1], axis=0) # Pone el fondo azul
	for i in range(fixedcolors+1,ncolors):
		gris = [240/255, 235/255, 235/255, 1]
		color_list = np.insert(color_list,i,gris, axis=0) # Pone el resto en gris
	return base.from_list('rainbow5', color_list, ncolors)


def plotheat(data, oname, nc, t):
	mask_ut=np.triu(np.ones(data.shape)).astype(np.bool)
	# ax = sns.heatmap(A,xticklabels=False, yticklabels=False, cmap="Accent_r")
	name = oname + "-heat.png"
	print(name)
	# base = plt.cm.get_cmap('rainbow') #paired
	# fixedcolors =10
	# color_list = base(np.linspace(0, 1, fixedcolors))
	# # color_list = np.insert(color_list,0,[1, 1, 1, 1], axis=0) # Pone el fondo blanco
	# color_list = np.insert(color_list,0,[21/255, 25/255, 235/255, 1], axis=0) # Pone el fondo azul

	# for i in range(fixedcolors+1,nc):
	# 	gris = [240/255, 235/255, 235/255, 1]
	# 	color_list = np.insert(color_list,i,gris, axis=0) # Pone el resto en gris
	# # print(color_list)

	# sns.heatmap(A, center=0, cmap=sns.diverging_palette(220, 20, as_cmap=True))	
	ax = sns.heatmap(data, mask=mask_ut, xticklabels=False, yticklabels=False, 
		cmap=cmapfixed(nc))
	plt.title(t)
	# plt.ylabel('Gene start position')
	# plt.xlabel('Gene start position')
	plt.tight_layout()
	plt.savefig(name,dpi=300)
	plt.clf()
	plt.close() # to clean memory
	# plt.show()


def plotcum(data, oname, nc, t):
	print("Plotting cum")
	pd.set_option('display.max_rows', None)
	name = oname + "-cum.png"
	print(data)
	data = data.reset_index(drop = True)
	cnames = data["clusterid"].unique()
	cums = pd.DataFrame(columns = cnames)
	for cid in cnames: # cid = [c14,c2,c4,...]
		cums[cid] = data['clusterid']
		cums.loc[cums[cid] != int(cid), cid] = 0  
		cums.loc[cums[cid] == int(cid), cid] = 1
		cums[cid] = cums[cid].cumsum()
	my_cmap = matplotlib.colors.ListedColormap(my_rgbs, name='my_colormap_name')
	cums.plot(figsize=(9,7), title = t)
	plt.show()
	print(cums)
	sys.exit(15)

		
		
		
		
		# cums[cid].plot() # Plot of every cluster
	# print(t.clusterid.value_counts())
	
	
	# plt.xlabel("Gene start")
	# plt.ylabel("Number of Genes")
	plt.legend(frameon=False, loc='upper left', ncol=3, fontsize = 'x-small')	
	oname = plotsdir + 'clustercum-chr' + chrid + ".png"
	# print("Writing plot: " + oname)
	# plt.savefig(oname,dpi=300)
	# plt.clf()
	# plt.close() # to clean memory
	plt.show()


def addcolorcol(df, colname="clusterid"):
	# Suma una columna de color con el value count de otra columna
	vc = df[colname].value_counts()
	vc = vc.rename_axis('unique_values').reset_index(name='counts')
	print(vc["unique_values"].tolist())
	# df = df.set_index('clusterid')
	# df = df.loc[vc["unique_values"].tolist()].reset_index() #order by list
	# print(df)
	df["color"] = 1
	df["color"] = df["color"].astype('int')
	i = 0
	for cl in vc["unique_values"].tolist():
		df.loc[ df[colname] == cl, "color"] = int(i)
		i = i + 1
	
	return df



############################################################
###  MAIN
############################################################

#! python py/03-corr-heatmaps.py Data/clustering/Healthy-clusters.tsv
# fname = sys.argv[1] # tsv filename: gen,clusterid,cost,genestart,chromid

dirpng = "Plots"
# healthyDir = "Data/Clustered/Healthy/Healthy-chr1-clusters.tsv"
Dir = "Data/Clustered/"

chrs = np.arange(1, 22).tolist()
chrs.append('X')
chrs = ['19'] 

for chrom in chrs:
	
	subtype = "Healthy"
	fname = Dir + "/" + subtype + "/" + subtype + "-chr" + str(chrom) + "-clusters.tsv"
	d = dirpng + "/chr" + str(chrom)
	logprint("Plot dir: %s" % d)
	logprint("Plot label: %s" % subtype)
	logprint("Plot chr name: %s" % chrom)
	fp = Path(d)
	fp.mkdir(parents=True, exist_ok=True)


	#! 1. PREPARING DATA
	df = pd.read_csv(fname, sep = ",")
	df = df.sort_values(by=['chromosome','gstart']) # Order by cromosome and clusterid
	nclusters = df["clusterid"].max()	
	print("Number of Clusters: %s" % nclusters)
	print("Number of Genes: %s" % df.shape[0]) # Print rows	

	df = addcolorcol(df)
	A = list2adj(df, start=5)

	# ! 2. PLOTTING
	oname = d + "/" + subtype + "-chr" + str(chrom) + "-gstart"
	print("Plot name: " + oname)
	title = 'Clusters in whole genome of ' + subtype + " chr " + str(chrom)
	plotheat(A, oname, nclusters, title)
	title = 'Cumulative Distribution Chromosome ' + str(chrom)
	# plotcum(df, oname, nclusters, title)
	sys.exit(15)
	
	
	subtypes = ["Basal","Her2","LumA","LumB"] 
	for subtype in subtypes:
		logprint("Subtype: %s" % subtype)
		fname = "Data/Clustered/" + subtype + "/" + subtype + "-chr" + str(chrom) + "-clusters.tsv"
		print(fname)
		df = pd.read_csv(fname, sep = ",")
		df = df.sort_values(by=['chromosome','gstart'])

		df = addcolorcol(df)
		A = list2adj(df, start=5)

		oname = d + "/" + subtype + "-chr" + str(chrom) + "-gstart" #+ "-clusters.png"
		# print(oname)
		title = 'Clusters in whole genome of ' + subtype + " chr " + str(chrom)
		plotheat(A, oname, nclusters, title)




