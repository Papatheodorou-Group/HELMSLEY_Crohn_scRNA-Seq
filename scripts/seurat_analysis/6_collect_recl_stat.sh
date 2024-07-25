
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "2" ]]
then
    echo ""
    echo "Usage: bash $0 <param_dir> <out_dir>"
    echo ""
    exit
fi

PARAM_DIR="$1"
OUT_DIR="$2"
GCA_STATS="$3"

PARAM="${PARAM_DIR}/param_clust_filt.txt"
PCA_DIR="${OUT_DIR}/hvg_pca_clust"
SILH_DIR="${OUT_DIR}/silhouette"
STAT_DIR="${OUT_DIR}/stats"
OBJ="${PCA_DIR}/object.Rds"
OUT="${OUT_DIR}/STATS.tsv"

Rscript "${SCRIPT_DIR}/collect_cl_stat.R" ${PARAM} "${PCA_DIR}" "${SILH_DIR}" "${OBJ}" "${OUT}"
# Rscript "${SCRIPT_DIR}/confusion_matrices.R" ${PARAM} "${STAT_DIR}" "${OBJ}"
Rscript "${SCRIPT_DIR}/confusion_matrices_v2.R" ${PARAM} "${STAT_DIR}" "${OBJ}" # NEW: change the metadata field according to the new way to add
                                                                                #      doublet information (see add_metadata_df.sh)
Rscript "${SCRIPT_DIR}/plot_confusion_matrices.R" ${PARAM} "${STAT_DIR}" "${GCA_STATS}"

exit


