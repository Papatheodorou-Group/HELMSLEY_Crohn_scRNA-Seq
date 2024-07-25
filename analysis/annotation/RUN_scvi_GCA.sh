SCRIPT_DIR="../../scripts/scvi"
SEURAT_SCRIPT_DIR="../../scripts/seurat_analysis"
DATA_DIR="../../data/GCA"

SEURAT_DIR="../gene_expression/merged"
MATRIX_DIR="../cell_calling"
Q_DATA_DIR="../../data/"

MEM="20000"
CORES="2"

############### GCA ################

OUT_MODEL_DIR="GCA/scvi_models"
OUT_AN_MODEL_DIR="GCA/scanvi_models"
OUT_OBJ_DIR="GCA/adata"
OUT_PLOT_DIR="GCA/plots"
OUT_STAT_DIR="GCA/stats"

IN_OBJ="${DATA_DIR}/obj_healthy_adult_pediatric.h5ad"

SET="adult_pediatric"

SCVI_REDUCT="X_scVI"
SCANVI_REDUCT="X_scANVI_category"

#####################################

[ -d "${OUT_MODEL_DIR}" ] || mkdir -p "${OUT_MODEL_DIR}"
[ -d "${OUT_OBJ_DIR}" ] || mkdir -p "${OUT_OBJ_DIR}"
[ -d "${OUT_STAT_DIR}" ] || mkdir -p "${OUT_STAT_DIR}"

# generate scVI and scANVI models on GCA

MODEL="${OUT_MODEL_DIR}/${SET}"
MODEL_AN="${OUT_AN_MODEL_DIR}/${SET}"
OUT_OBJ="${OUT_OBJ_DIR}/${SET}.h5ad"
STAT_R_PREFIX="${OUT_STAT_DIR}/${SET}"

OUT_DIR="${OUT_PLOT_DIR}/${SET}"
[ -d "${OUT_DIR}" ] || mkdir -p "${OUT_DIR}"

R_OBJ="${OUT_OBJ}"

# train the models and save them on anndata structure
ID1="step1"
COMMAND="${SCRIPT_DIR}/create_models_on_ref.sh"
bsub -J ${ID1} -e ${ID1}.STDERR -o ${ID1}.STDOUT -M ${MEM} -n ${CORES} bash ${COMMAND} "${IN_OBJ}" "${MODEL}" "${MODEL_AN}" "${R_OBJ}" ${SCVI_REDUCT} ${SCANVI_REDUCT} "${STAT_R_PREFIX}"

# plot the cells on PCA or scvi latent space
ID2="step2"
COMMAND="${SCRIPT_DIR}/plot_celltypes.sh"
bsub -J ${ID2} -w "done(${ID1})" -e ${ID2}.STDERR -o ${ID2}.STDOUT -M ${MEM} bash ${COMMAND} "${R_OBJ}" "${OUT_DIR}" ${SCVI_REDUCT}

# plot the cells on PCA or scanvi latent space
ID3="step3"
COMMAND="${SCRIPT_DIR}/plot_celltypes.sh"
bsub -J ${ID3} -w "done(${ID1})" -e ${ID3}.STDERR -o ${ID3}.STDOUT -M ${MEM} bash ${COMMAND} "${R_OBJ}" "${OUT_DIR}" ${SCANVI_REDUCT}

exit

