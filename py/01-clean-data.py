import sys, os
import pandas as pd
import shutil

# fname = "Data/Basal-All.txt"
fname = sys.argv[1]
print("Expression Matrix: ", fname)
cleandir = "Data/clean"

# File mods
if not os.path.exists(cleandir):
	os.mkdir(cleandir)
	print("Creating [" + cleandir + "]")

l = fname.split('/')[2].split('-')[0]
print(l)
oname = cleandir + "/" + l + "-clean.tsv" 

# 1) READ EXPRESSION DATA
matexp = pd.read_csv(fname, sep = "\t", header=None)
matexp = matexp.rename(columns={ matexp.columns[0]: "gname" })
print("Original Expression Matrix: ",matexp.shape)
#? EXAMPLE OF BAD ROW: gname = 3974.44679218112
matexp.dropna(inplace = True) # CLEAN NaN IN BAD ROWS 
# EXAMPLE OF REGEX IN COLUMN
# patt = matexp['gname'].str.contains('^[1-9]+')
# pd.set_option('display.max_rows', None) ## PRINT ALL ROWS IN PANDAS
# print(matexp[patt])
matexp.reset_index(drop=True,inplace=True)
# gname_col = matexp.iloc[:,[0]] ## CUT GENE NAMES
# matexp = matexp.iloc[:, 1:] ## CUT ONLY VALUES
matexp = matexp.iloc[:, 0:-1] ## CUT ERROR COLUMN
print("Clean Expression Matrix: ",matexp.shape)
# print("Clean Genes: ",gname_col.shape)

# 2) READ BIOMART DATA
mart = pd.read_csv("Data/biomart-20190718.txt", sep = "\t")
mart.drop(['Gene stable ID','Gene stable ID version'], inplace=True, axis = 1)
# Clean only standard chromosomes
patt = mart['Chromosome/scaffold name'].str.contains('^[1-9]+|^[XY]') 
mart = mart[patt]
#? Genes with same chrom & name but differente gstart. Ex,PRSS50
# Clean genes with same chrom & name but differente gstart
mart.drop_duplicates(subset=['Gene name'], inplace=True)  

# 3) MERGE BIOMART & EXPRESSION DATA
matexp = pd.merge(matexp, mart, left_on='gname', right_on='Gene name', how='left')
matexp.drop(['Gene name'], inplace=True, axis = 1)
matexp.rename(columns={'Gene start (bp)': 'gstart', 'Chromosome/scaffold name': 'chromname'}, inplace=True)
print("Merge: ", matexp.shape)
## Genes in Expression not in Biomart
print("Genes without Biomart: ", matexp[matexp.isna().any(axis=1)].shape[0])
matexp.dropna(inplace = True) # Clean NaNs en Biomart genes
matexp = matexp.astype({"gstart": int})
matexp.reset_index(drop=True,inplace=True)
print("Clean matexp", matexp.shape)
matexp.sort_values(['chromname', 'gstart'])
matexp.to_csv(oname,index=False,sep='\t')

# print("Saving by chromosomes ...")
# chroms = matexp["chromname"].value_counts().index
# for c in chroms:
# 	matbychrom = matexp[matexp['chromname'] == c]
# 	oname = cleandir  + "/" + l
# 	oname = oname  + "-clean-chr" + c + ".tsv" 	
# 	matbychrom.to_csv(oname,index=False,header=False,sep='\t')

