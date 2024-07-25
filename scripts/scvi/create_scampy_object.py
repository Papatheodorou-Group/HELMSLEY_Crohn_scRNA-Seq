# source ~/.bashrc
# conda activate sklearn-env

from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <matrix_dir> <metadata> <feature_mapping> <out>\n"
arg1 = "\nmatrix_dir: directory containing matrix.mtx, genes.tsv, barcodes.tsv"
arg2 = "\nmetadata: tsv file containing additional information"
arg3 = "\nfeature_mapping: tsv file containing gene symbols (with Seurat suffix) and ENSEMBL id"
arg4 = "\nout: file containing scanpy object"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3 + arg4)
(options, args) = parser.parse_args()

if len(args) < 4:
	sys.stderr.write("ERROR: <matrix_dir>, <metadata>, <feature_mapping> and <out> are required\nTry --help for help\n")
	exit()
	
MATRIX_DIR=args[0]
META=args[1]
FEATURE=args[2]
OUT=args[3]

###################
# PROCESSING
###################

import scanpy as sc
import pandas as pd
import numpy as np
from scipy import sparse
import pickle

# load the cell-gene matrix and the annotations
M = sc.read(MATRIX_DIR + "/matrix.mtx") # or scanpy.read_mtx
genes = pd.read_csv(MATRIX_DIR + "/genes.tsv", header = None)
barcodes = pd.read_csv(MATRIX_DIR + "/barcodes.tsv", header = None)

M.var = genes
ncells = M.shape[0] 
ngenes = M.shape[1]

M.obs = pd.read_csv(META, sep = "\t")
M.obs['Sample name'] = M.obs['Sample name'].astype(str)
M.obs['Diagnosis'] = M.obs['Diagnosis'].astype(str)
M.obs['Region code'] = M.obs['Region code'].astype(str)
		
feature = pd.read_csv(FEATURE, sep = "\t")
		
# match features with genes in the matrix
f_bool = np.isin(feature["symbol"], genes)
f_idx = np.where(f_bool)[0] # genes that are also in the matrix
idx = [] # index of common genes in the matrix
for i in range(feature.shape[0]):
	if feature["symbol"][i] in genes.values:
		c_idx = genes[genes[0] == feature["symbol"][i]].index[0]
		idx.append(c_idx)

data_dict = {"gene_ids": ['nan'] * ngenes, "symbol": genes[0].values, "feature_types": ['Gene expression'] * ngenes} 
M.var = pd.DataFrame.from_dict(data_dict)
M.var.iloc[idx,:] = feature.iloc[f_idx,:]
M.var.set_index('symbol', inplace=True)

M.layers["counts"] = M.X.copy()  # preserve counts
sc.pp.normalize_total(M, target_sum=1e6)
sc.pp.log1p(M)
M.raw = M

M.write_h5ad(OUT)


exit()

