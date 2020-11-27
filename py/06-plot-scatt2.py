import pandas as pd
import numpy as np
import seaborn as sns
import glob, sys, os 
import matplotlib.pyplot as plt
from termcolor import cprint

logprint = lambda x: cprint(x, 'red', attrs=["bold"])

corrsdir = "Data/Pearson/"

for f in glob.glob(corrsdir + '*.csv'): 

	subtype = f.split('/')[-1].split('-')[0]
	logprint(subtype)
	print("Loding data ...")
	df = pd.read_csv(f, sep = ",")
	df['src_chrom']= df['src_chrom'].astype(str, copy=False)
	df['dst_chrom']= df['dst_chrom'].astype(str, copy=False)
	df.drop(['src', 'dst'], axis=1, inplace=True)

	print("Plot chromosomes ...")
	
	for chrom in df.src_chrom.unique():
		logprint(chrom)
		corr = df[df['src_chrom']==chrom]
		corr = corr[corr['dst_chrom']==chrom]
		corr = corr.assign(distance = corr['dst_gstart'] - corr['src_gstart'])
		corr.sort_values(by=['distance'], inplace=True)
		# print(corr)
		
		corr.plot.scatter(x='distance', y='pearson',
			title= "Pearson Distance", c='Black', s=2) 
		plt.axhline(y=0, color='r', linestyle='-', markersize=20)

		fout = 'Plots/chr' + chrom + '/'
		os.makedirs(fout, exist_ok=True)
		fout =  fout + subtype + '-chr' + chrom + '-pdist.png'
		print(fout)

		plt.savefig(fout,dpi=300)
		plt.clf()
		plt.close()

		# input("Press the <ENTER> key to continue...")
		# sys.exit(15)









	# for chrom, g in df.groupby(['src_chrom']):
	# 	logprint(chrom)
	# 	corr = g[g['dst_chrom']==chrom]
	# 	corr = corr.assign(distance = corr['dst_gstart'] - corr['src_gstart'])
	# 	corr.sort_values(by=['distance'], inplace=True)
	# 	# print(corr)
	# 	corr.plot.scatter(x='distance', y='pearson', 
	# 		title= "Pearson Distance", c='Black', s=2);
	# 	plt.axhline(y=0, color='r', linestyle='-', markersize=20)
	# 	# plt.show()
	# 	fout = 'Plots/chr' + chrom + '/'
	# 	os.makedirs(fout, exist_ok=True)
	# 	fout =  fout + subtype + '-chr' + chrom + '-pdist.png'
	# 	print(fout)
	# 	plt.savefig(fout,dpi=300)
	# 	plt.clf()
	# 	plt.close() # to clean memory

	# 	logprint("Delete corr")
	# 	del corr

	# logprint("Delete df")
	# del df

#64 + 78