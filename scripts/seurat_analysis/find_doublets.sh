
WD=$( cd $( dirname $0 ) ; pwd )

if [[ $# -lt "6" ]]
then
    echo ""
    echo "Usage: bash $0 <script_dir> <IN_DIR> <OUT_DIR> <param_dir> <doublet_perc> <slot>"
    echo ""
    exit
fi

SCRIPT_DIR="$1"
IN_DIR="$2"
OUT_DIR="$3"
PARAM_DIR="$4"
DOUBLET_PERC="$5"
SLOT="$6"

R_SCRIPT_HVG_PCA="${SCRIPT_DIR}/Seurat_hvg_pca.R"
R_SCRIPT_DOUBLETS="${WD}/find_doublets.R"
R_SCRIPT_DOUBLETS_STATS="${WD}/collect_doublet_stats.R"

HPC_DIR="${OUT_DIR}/hvg_pca_clust"
PLOT_DIR="${OUT_DIR}/plots"

IN_OBJ="${OUT_DIR}/object.Rds"
OUT_OBJ="${HPC_DIR}/object.Rds"

ORIG_DIR="${IN_DIR}/hvg_pca_clust"
ORIG_OBJ="${ORIG_DIR}/object.Rds"

STAT_PREFIX="${OUT_DIR}/doublet"

PARAMS_PCA="${PARAM_DIR}/param_hvg_pca_filt.txt"

# run feature selection and PCA
# do not run Jackstraw test here!
echo "run feature selection and PCA"
Rscript "${R_SCRIPT_HVG_PCA}" \
        "${IN_OBJ}" \
        "${HPC_DIR}" \
        "${PARAMS_PCA}" \
        1> "${OUT_DIR}/hvg_pca.STDOUT" \
        2> "${OUT_DIR}/hvg_pca.STDERR"

# predict doublets
echo "predict doublets"
Rscript "${R_SCRIPT_DOUBLETS}" \
        "${OUT_OBJ}" \
        "${HPC_DIR}" \
        "${PLOT_DIR}" \
        "${PARAMS_PCA}" \
        ${DOUBLET_PERC} \
        1> "${OUT_DIR}/doublets.STDOUT" \
        2> "${OUT_DIR}/doublets.STDERR"

echo "collect doublets stats"
Rscript "${R_SCRIPT_DOUBLETS_STATS}" \
        "${ORIG_OBJ}" \
        "${SLOT}" \
        "${OUT_OBJ}" \
        "${STAT_PREFIX}" \
        1> "${OUT_DIR}/doublets_stats.STDOUT" \
        2> "${OUT_DIR}/doublets_stats.STDERR"


exit
