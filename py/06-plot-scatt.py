import pandas as pd
import numpy as np
import glob, sys, os
import matplotlib.pyplot as plt
from termcolor import cprint
from pandarallel import pandarallel

pandarallel.initialize()

def randrow(row):
	np.random.shuffle(row)
	return row

def getcorr(exprmat):
	corr = exprmat.iloc[:,1:exprmat.shape[1]-2]
	corr = corr.T.corr()
	# corr['gname'] = exprmat['gname']
	corr.index = list(exprmat["gname"])
	corr.columns = list(exprmat["gname"])
	m,n = corr.shape
	corr[:] = np.where(np.arange(m)[:,None] >= np.arange(n),np.nan,corr)
	corr = corr.stack().reset_index()
	corr.columns = ["src","dst","cor"]
	corr = pd.merge(corr, gstart, how='left', left_on='src',right_on='gname')
	corr.drop('gname', axis = 1, inplace = True)
	corr.rename(columns = {'gstart':'src_gstart'}, inplace = True)
	corr = pd.merge(corr, gstart, how='left', left_on='dst',right_on='gname')
	corr.drop('gname', axis = 1, inplace = True)
	corr.rename(columns = {'gstart':'dst_gstart'}, inplace = True)
	corr = corr.assign(distance = corr['dst_gstart'] - corr['src_gstart'])
	corr.sort_values(by=['distance'], inplace=True)
	print(corr)
	return corr

logprint = lambda x: cprint(x, 'red', attrs=["bold"])

corrsdir = "Data/clean/"

for f in glob.glob(corrsdir + '*.tsv'):

	subtype = f.split('/')[-1].split('-')[0]
	logprint(subtype)
	print("Loding data ...")
	df = pd.read_csv(f, sep = "\t")
	df['chromname']= df['chromname'].astype(str, copy=False)

	print("Processing Chromosomes ...")
	for chrom in df.chromname.unique():
		logprint("Chromosome: " + chrom)
		expr = df[df['chromname']==chrom]
		gstart = expr.loc[:,['gname','gstart']] # Table:gname,gstart

		# Calculation Pearson standard
		corr1 = getcorr(expr)
		ax1 = corr1.plot.scatter(x='distance', y='cor', s=4, linewidths=0.2,
			color='r', edgecolors='#000000')
		
		# Calculation Pearson with Shuffle by cols 
		expr2 = expr.iloc[:,1:expr.shape[1]-2]#.to_numpy()
		expr2 = expr2.parallel_apply(randrow, axis=1)
		expr2 = expr2.join(expr.gname)
		expr2 = expr2.join(expr.gstart)
		corr2 = getcorr(expr2)
		corr2.plot.scatter(x='distance', y='cor', s=4, linewidths=0.2, 
			color='b', edgecolors='#000000', ax=ax1)

		# Save Plot
		ax1.set_ylim(-1, 1)		
		fout = 'Plots/chr' + chrom + '/'
		os.makedirs(fout, exist_ok=True)
		fout =  fout + subtype + '-chr' + chrom + '-pdist.png'
		print(fout)
		plt.savefig(fout,dpi=300)
		plt.clf() 
		plt.close()
		# plt.show()

	





























	# print(df)
	# expr = df.iloc[:,1:df.shape[1]-2]
	# print(expr)
	# # 
	
	# expr = expr.sample(frac=1, axis=1).sample(frac=1).reset_index(drop=True)
	# print(expr)
	
	# break

	# print("Plot chromosomes ...")
	
	# for chrom in df.src_chrom.unique():
	# 	logprint(chrom)
	# 	corr = df[df['src_chrom']==chrom]
	# 	corr = corr[corr['dst_chrom']==chrom]
	# 	corr = corr.assign(distance = corr['dst_gstart'] - corr['src_gstart'])
	# 	corr.sort_values(by=['distance'], inplace=True)
	# 	# print(corr)
		
	# 	corr.plot.scatter(x='distance', y='pearson',
	# 		title= "Pearson Distance", c='Black', s=2) 
	# 	plt.axhline(y=0, color='r', linestyle='-', markersize=20)

	# 	fout = 'Plots/chr' + chrom + '/'
	# 	os.makedirs(fout, exist_ok=True)
	# 	fout =  fout + subtype + '-chr' + chrom + '-pdist.png'
	# 	print(fout)

	# 	plt.savefig(fout,dpi=300)
	# 	plt.clf()
	# 	plt.close()

	# 	# input("Press the <ENTER> key to continue...")
	# 	# sys.exit(15)









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