
from optparse import OptionParser
import os
import sys

MIN_CELLS = 5

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <in_obj1> <in_obj2> <outdir> <reduct_name> <meta>\n"
arg1 = "\nin_obj1: input anndata object (query)"
arg2 = "\nin_obj2: input anndata object (reference)"
arg3 = "\nout_dir: output directory where to save the plots"
arg4 = "\nreduct_name: name of the dimensional reduction to use"
arg5 = "\nbroad meta: for plotting and subsetting"
arg6 = "\nfine meta: for plotting"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3 + arg4 + arg5 + arg6)
(options, args) = parser.parse_args()

if len(args) < 4:
	sys.stderr.write("ERROR: <in_obj1>, <in_obj2>, <out_dir>, <reduct_name> and <meta> are required\nTry --help for help\n")
	exit()
	
Q_OBJ=args[0]
R_OBJ=args[1]
DIR=args[2]
REDUCT_NAME=args[3]
BROAD_META=args[4]
FINE_META=args[5]

###################
# EXECUTE
###################

import scanpy as sc
from matplotlib import pyplot as plt
import numpy as np

adata_query = sc.read_h5ad(Q_OBJ)
adata_ref = sc.read_h5ad(R_OBJ)
adata_full = adata_query.concatenate(adata_ref)

adata_full.obs['batch'] = adata_full.obs.batch.cat.rename_categories(["Query", "Reference"])

sc.pp.neighbors(adata_full, use_rep=REDUCT_NAME)
sc.tl.umap(adata_full, min_dist=0.3)

sc.pl.umap(
    adata_full,
    color=[BROAD_META, "batch"],
    frameon=False,
    show=False
)
plt.savefig(DIR + "/UMAP_all_" + REDUCT_NAME + "_" + BROAD_META + ".pdf", format="pdf", bbox_inches='tight')

ax = sc.pl.umap(
    adata_full,
    frameon=False,
    show=False
)
sc.pl.umap(
    adata_full[: adata_query.n_obs],
    color=["Sample name"],
    frameon=False,
    ax=ax
)
plt.savefig(DIR + "/UMAP_query_" + REDUCT_NAME + "_" + BROAD_META + ".pdf", format="pdf", bbox_inches='tight')
plt.close()

categories = np.unique(adata_full.obs[BROAD_META])
for c in categories:

	if sum(adata_full.obs[BROAD_META] == c) > MIN_CELLS:

		ax = sc.pl.umap(
		    adata_full,
		    frameon=False,
		    show=False
		)
		sc.pl.umap(
		    adata_full[adata_full.obs[BROAD_META] == c],
		    color=[FINE_META],
		    frameon=False,
		    ax = ax
		    )
		plt.savefig(DIR + "/UMAP_all_" + REDUCT_NAME + "_" + c + ".pdf", format="pdf", bbox_inches='tight')
		plt.close()
	
		ax = sc.pl.umap(
		    adata_full,
		    frameon=False,
		    show=False
		    )
		sc.pl.umap(
	    	adata_full[(adata_full.obs[BROAD_META] == c) & (adata_full.obs["batch"] == "Query")],
	    	    color=[FINE_META],
	    	    frameon=False,
	    	    ax=ax
		)
		plt.savefig(DIR + "/UMAP_query_" + REDUCT_NAME + "_" + c + ".pdf", format="pdf", bbox_inches='tight')
		plt.close()

exit()

