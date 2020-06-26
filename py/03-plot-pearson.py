import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from palettable.colorbrewer.diverging import RdBu_5
from termcolor import colored, cprint
import sys
logprint = lambda x: cprint(x, 'red', attrs=["bold"])


#############################
# MAIN
#############################
print("Processing input")
fname = sys.argv[1] # tsv filename: gen,clusterid,cost,genestart,chromid
dirpng = "Plots/"
label = fname.split('/')[-1].split('-')[0]
dirpng = dirpng + label
path = Path(dirpng).mkdir(parents=True, exist_ok=True)
logprint("Plot dir: %s" % dirpng)
logprint("Label: %s" % label)

# ALL CHROMOSOME WORK
print("Getting correlation")
df = pd.read_csv(fname, sep = "\t")
cor =  df.iloc[:,1:df.shape[1]-2].T.corr()

print("Ploting")
plt.figure(figsize=(7,7),dpi=300)
masking = np.tril(cor)
sns.heatmap(cor ,xticklabels=False, yticklabels=False, 
	mask=masking, vmin=-1., vmax=1., 
	cmap=RdBu_5.mpl_colormap)
plt.title('Pearson correlation' + label)
plt.ylabel('Gene start position')
plt.xlabel('Gene start position')
plt.tight_layout()
pngname = dirpng + '/' + label + '-allchroms.png'
plt.savefig(pngname)
plt.clf()
plt.close() # to clean memory


# BY CROMOSOME WORK
df['chromname']= df['chromname'].astype(str)
for gr_name, df_chr in df.groupby('chromname'):
	print(gr_name, len(df_chr.index))
	print("Working in chr: ", gr_name)
	d = df_chr.sort_values('gstart')
	ncols = df.shape[1]-2
	cor = d.iloc[:,1:ncols].T.corr()
	print("Ploting chr: ", gr_name)
	plt.figure(figsize=(5,5),dpi=300)
	masking = np.tril(cor)
	sns.heatmap(cor ,xticklabels=False, yticklabels=False, 
		mask=masking, vmin=-1., vmax=1., 
		cmap=RdBu_5.mpl_colormap)
	t = "Pearson correlation: " + label 
	t = t + ', chr: ' + str(gr_name)
	t = t + ', genes: ' + str(len(d.index))
	plt.title(t)
	plt.ylabel('Gene start position')
	plt.xlabel('Gene start position')
	plt.tight_layout()
	pngname = dirpng + '/' + label
	pngname = pngname + '-chr' + str(gr_name) + '.png'
	plt.savefig(pngname)
	plt.clf()
	plt.close() # to clean memory


