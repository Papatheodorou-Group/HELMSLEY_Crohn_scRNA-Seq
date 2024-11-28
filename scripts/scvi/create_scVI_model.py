from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

BATCH_KEY = "Diagnosis" # for GCA: put this instead of "batch", because apparently there are too few cells for some batches, which end up in errors (ValueError: b'There are other near singularities as well. 0.090619\n')

usage = "python %prog [options] <in_obj> <out> <out_obj> <reduct_name>\n"
arg1 = "\nin_obj: input h5ad object"
arg2 = "\nout: output scvi model"
arg3 = "\nout_obj: output h5ad object"
arg4 = "\nreduct_name: name of the dimensional reduction (e.g. scvi, scvi_TIL...)"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3 + arg4)
(options, args) = parser.parse_args()

if len(args) < 4:
	sys.stderr.write("ERROR: <in_obj>, <out>, <out_obj> and <reduct_name> are required\nTry --help for help\n")
	exit()
	
IN_OBJ=args[0]
OUT=args[1]
OUT_OBJ=args[2]
REDUCT_NAME=args[3]
if len(args) > 4:
	BATCH_KEY=args[4]

###################
# EXECUTE
###################

import scanpy as sc
import scvi

# UserWarning: Since v1.0.0, scvi-tools no longer uses a random seed by default
scvi.settings.seed = 0 

adata = sc.read_h5ad(IN_OBJ)

adata.layers["counts"] = adata.X.copy()  # preserve counts
sc.pp.normalize_total(adata, target_sum=1e6)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(
    adata,
    n_top_genes=5000, 
    subset=True,
    layer="counts",
    flavor="seurat_v3",
    batch_key=BATCH_KEY
)

# NotImplementedError: scArches currently does not support models with extra categorical covariates.
# -> hence leave "categorical_covariate_keys" empty 
scvi.model.SCVI.setup_anndata(
    adata,
    layer="counts",
#    categorical_covariate_keys=["Sample name", "Diagnosis", "Region code"],
#    continuous_covariate_keys=["pct_counts_mt"] # remove this so I do not have to compute it for the query
)

arches_params = dict(
    use_layer_norm="both",
    use_batch_norm="none",
    encode_covariates=True,
    dropout_rate=0.2,
    n_layers=2,
)

model = scvi.model.SCVI(adata, **arches_params)
model.train()
# model.train(max_epochs = 2) # FOR DEBUGGING

model.save(OUT, overwrite = True)

adata.obsm["X_" + REDUCT_NAME] = model.get_latent_representation()
denoised = model.get_normalized_expression(adata, library_size=1e6)
adata.layers["normalized_" + REDUCT_NAME] = model.get_normalized_expression(library_size=10e6)

adata.write_h5ad(OUT_OBJ)

exit()

