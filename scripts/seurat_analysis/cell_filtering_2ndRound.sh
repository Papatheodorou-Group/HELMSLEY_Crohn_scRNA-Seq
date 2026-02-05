
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt 3 ]]
then
    echo "Usage: bash $0 <dir2> <dr> <outdir> <dir1> <mode> <cl>"
    exit
fi

dir_doublets=$1 # dir where doublet labelling is stored
dr=$2 
outdir=$3
if [[ $# -gt 4 ]]
then
    dir_clusters=$4 # dir where clustering is stored
    mode=$5
    cl=$6
fi

[[ -d ${outdir} ]] || mkdir -p ${outdir}
cells="${outdir}/cells.txt"

# extract the cells
if [[ $# -lt 6 ]]
then
    Rscript ${SCRIPT_DIR}/cell_subset_doublet.R ${dir_doublets}/hvg_pca_clust/object.Rds ${dr} ${cells} 
else
    Rscript ${SCRIPT_DIR}/cell_subset_doublet.R ${dir_doublets}/hvg_pca_clust/object.Rds ${dr} ${cells} ${dir_clusters}/hvg_pca_clust/object.Rds ${mode} ${cl}
fi

# subset the object 
Rscript ${SCRIPT_DIR}/obj_subset.R ${dir_doublets}/hvg_pca_clust/object.Rds ${cells} ${outdir}/object.Rds

exit

