import pandas as pd
import numpy as np
import seaborn as sns
import glob, sys
import matplotlib.pyplot as plt

logprint = lambda x: cprint(x, 'red', attrs=["bold"])

corrsdir = "Data/Pearson/"

for f in glob.glob(corrsdir + '*.csv'): 

	subtype = f.split('/')[-1].split('-')[0]
	print(subtype)	
	df = pd.read_csv(f, sep = ",")
	df['src_chrom']= df['src_chrom'].astype(str)
	df['dst_chrom']= df['dst_chrom'].astype(str)

	for chrom, g in df.groupby(['src_chrom']):
		print(chrom)
		corr = g[g['dst_chrom']==chrom]
		corr = corr.assign(distance = corr['dst_gstart'] - corr['src_gstart'])
		corr.sort_values(by=['distance'], inplace=True)
		print(corr)
		corr.plot.scatter(x='distance', y='pearson', 
			title= "Pearson Distance", c='Black');
		plt.axhline(y=0, color='r', linestyle='-', markersize=12)
		plt.show()
		fout = 'Plots/chr' + chrom + '/'
		fout =  fout + subtype + '-chr' + chrom + '-pdist.png'
		print(fout)
		# plt.savefig(fout,dpi=300)
		# plt.clf()
		# plt.close() # to clean memory