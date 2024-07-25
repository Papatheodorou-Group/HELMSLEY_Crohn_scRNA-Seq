
from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <in_obj> <dir>\n"
arg1 = "\nin_obj: input anndata object"
arg2 = "\ndir: output directory where to save the plots"
arg3 = "\nreduct_name: name of the dimensional reduction to use"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3)
(options, args) = parser.parse_args()

if len(args) < 3:
	sys.stderr.write("ERROR: <in_obj>, <dir> and <reduct_name> are required\nTry --help for help\n")
	exit()
	
IN_OBJ=args[0]
DIR=args[1]
REDUCT_NAME=args[2]

###################
# EXECUTE
###################

import scanpy as sc
from matplotlib import pyplot as plt

adata = sc.read_h5ad(IN_OBJ)

# run PCA then generate UMAP plots
sc.tl.pca(adata)
sc.pp.neighbors(adata, n_pcs=30, n_neighbors=20)
sc.tl.umap(adata, min_dist=0.3)

plt = sc.pl.umap(
    adata,
    color=["category","Integrated_05"],
    ncols=2,
    frameon=False,
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig(DIR + "/UMAP_PCA_celltypes.pdf", format="pdf")

plt = sc.pl.umap(
    adata,
    color=["Diagnosis", "Sample name", "Region code"],
    ncols=4,
    frameon=False, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig(DIR + "/UMAP_PCA_metadata.pdf", format="pdf")

# use scVI latent space for UMAP generation
sc.pp.neighbors(adata, use_rep=REDUCT_NAME)
sc.tl.umap(adata, min_dist=0.3)

plt = sc.pl.umap(
    adata,
    color=["category","Integrated_05"],
    ncols=2,
    frameon=False, 
    return_fig=True,
    show=False
)
plt.set_size_inches(8, 4)
plt.savefig(DIR + "/UMAP_" + REDUCT_NAME + "_celltypes.pdf", format="pdf")

plt = sc.pl.umap(
    adata,
    color=["Diagnosis", "Sample name", "Region code"],
    ncols=4,
    frameon=False, 
    return_fig=True,
    show=False
)
plt.set_size_inches(12, 4)
plt.savefig(DIR + "/UMAP_" + REDUCT_NAME + "_metadata.pdf", format="pdf")

exit()

