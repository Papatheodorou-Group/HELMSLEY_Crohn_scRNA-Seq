
import scanpy as sc

gca_in = "GCA/Full_obj_raw_counts_nosoupx_v2.h5ad"
gca_out1 = "GCA/obj_healthy_adult_pediatric.h5ad"
gca_out2 = "GCA/obj_healthy_adult_pediatric_TIL.h5ad"
gca_out3 = "GCA/obj_healthy_adult_pediatric_SCL.h5ad"

diagnoses_to_extract = ['Healthy adult','Pediatric healthy']

# see data/Gut_cell_atlas/intestinal_tract_ID_info.tsv
# regions_to_extract = ['DUO', 'DCL', 'TIL', 'ACL', 'APD', 'CAE', 'TCL', 'ILE2', 'SCL', 'ILE1', 'ILE', 'JEJ', 'REC', 'MLN'] 

adata = sc.read_h5ad(gca_in)
adata = adata[adata.obs['Diagnosis'].isin(diagnoses_to_extract)]
adata.write_h5ad(gca_out1)

regions_to_extract = ['TIL']
adata_filt = adata[adata.obs['Region code'].isin(regions_to_extract)]
adata_filt.write_h5ad(gca_out2)

regions_to_extract = ['SCL']
adata_filt = adata[adata.obs['Region code'].isin(regions_to_extract)]
adata_filt.write_h5ad(gca_out3)

exit()


