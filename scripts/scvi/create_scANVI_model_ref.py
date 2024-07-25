from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <r_obj> <r_mod_in> <r_mod_out> <meta_for_scanvi> <reduct_name>\n"
arg1 = "\nr_obj_prefix: input h5ad object (without extension; reference)"
arg2 = "\nr_mod_in: input scvi model (reference)"
arg3 = "\nr_mod_out: output scanvi model (reference)"
arg4 = "\nmeta_for_scanvi: name of the metadata field to be used for training scANVI (e.g. category, cell type...)"
arg5 = "\nreduct_name: reduction name (latent space) (e.g. scANVI, scANVI_TIL...)"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3 + arg4 + arg5)
(options, args) = parser.parse_args()

if len(args) < 5:
	sys.stderr.write("5RROR: <r_obj_prefix>, <r_mod_in>, <r_mod_out>, <meta_for_scanvi> and <reduct_name> are required\nTry --help for help\n")
	exit()
	
R_OBJ_PREFIX=args[0]
R_MOD_IN=args[1]
R_MOD_OUT=args[2]
META=args[3]
REDUCT_NAME=args[4]

###################
# EXECUTE
###################

import scanpy as sc
import scvi

R_OBJ = R_OBJ_PREFIX + ".h5ad"
R_LAB = R_OBJ_PREFIX + "_" + META + ".tsv"
LATENT_ID = "X_" + REDUCT_NAME + "_" + META

# UserWarning: Since v1.0.0, scvi-tools no longer uses a random seed by default
scvi.settings.seed = 0 

adata_ref = sc.read_h5ad(R_OBJ)
vae_ref = scvi.model.SCVI.load(R_MOD_IN, adata=adata_ref)

adata_ref.obs["labels_scanvi"] = adata_ref.obs[META].values

vae_ref_scan = scvi.model.SCANVI.from_scvi_model(
    vae_ref, # tried with R_MOD_IN, as in the tutorial, but it doesn't work
    unlabeled_category="Unknown",
    labels_key="labels_scanvi"
)
vae_ref_scan.train(max_epochs=20, n_samples_per_label=100)
# vae_ref_scan.train(max_epochs=2, n_samples_per_label=100) # FOR DEBUGGING

vae_ref_scan.save(R_MOD_OUT, overwrite = True)

adata_ref.obsm[LATENT_ID] = vae_ref.get_latent_representation()
adata_ref.write_h5ad(R_OBJ)

labels = adata_ref.obs["labels_scanvi"]
labels.to_csv(R_LAB, sep = "\t") # NEW: print the labels (ground truth)

exit()

