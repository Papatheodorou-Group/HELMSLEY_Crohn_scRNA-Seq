from optparse import OptionParser
import os
import sys

###################
# PARSE OPTIONS
###################

usage = "python %prog [options] <q_obj> <r_mod_out> <meta_for_scanvi> <reduct_name>\n"
arg1 = "\nq_obj_prefix: input h5ad object file prefix (query)"
arg2 = "\nr_mod_out: scanvi model (reference)"
arg3 = "\nmeta_for_scanvi: name of the metadata field to be used for training scANVI (e.g. category, cell type...)"
arg4 = "\nreduct_name: reduction name (latent space) (e.g. scANVI, scANVI_TIL...)"
print("")

parser = OptionParser(usage=usage + arg1 + arg2 + arg3 + arg4)
(options, args) = parser.parse_args()

if len(args) < 4:
	sys.stderr.write("ERROR: <q_obj_prefix>, <r_mod_out>, <meta_for_scanvi> and <reduct_name> are required\nTry --help for help\n")
	exit()
	
Q_OBJ_PREFIX=args[0]
R_MOD_OUT=args[1]
META=args[2]
REDUCT_NAME=args[3]

###################
# EXECUTE
###################

import scanpy as sc
import scvi

Q_OBJ = Q_OBJ_PREFIX + ".h5ad"
SOFT_Q_LAB = Q_OBJ_PREFIX + "_soft" + META + ".tsv"
LATENT_ID = "X_" + REDUCT_NAME + "_" + META

# UserWarning: Since v1.0.0, scvi-tools no longer uses a random seed by default
scvi.settings.seed = 0 

adata_query = sc.read_h5ad(Q_OBJ)
lab_exists = "labels_scanvi" in adata_query.obs
if lab_exists:
    labels = adata_query.obs["labels_scanvi"] # NEW
adata_query.obs["labels_scanvi"] = "Unknown" # NEW

scvi.model.SCANVI.prepare_query_anndata(
    adata_query, 
    R_MOD_OUT
) 
vae_q = scvi.model.SCANVI.load_query_data(
    adata_query,
    R_MOD_OUT
)

# from scVI tutorial: "For training the query data, we recommend using a weight_decay of 0.0. This ensures the latent representation of the reference cells will remain exactly the same if passing them through this new query model."
vae_q.train(
    max_epochs=100,
#    max_epochs=10, # FOR DEBUGGING
    plan_kwargs={"weight_decay": 0.0}, 
    check_val_every_n_epoch=10
)

labels_hard = vae_q.predict()
labels_soft = vae_q.predict(soft = True)

adata_query.obsm[LATENT_ID] = vae_q.get_latent_representation()
adata_query.obs[META] = labels_hard

if lab_exists:
    adata_query.obs["labels_scanvi"] = labels # NEW
adata_query.write_h5ad(Q_OBJ)
labels_soft.to_csv(SOFT_Q_LAB, sep = "\t")

exit()

