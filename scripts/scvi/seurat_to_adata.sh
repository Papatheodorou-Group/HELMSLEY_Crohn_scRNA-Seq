source ~/.bashrc
conda activate scvi-env

# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

if [[ $# -lt "6" ]]
then
    echo "usage: bash seurat_to_adata.sh <seurat_script_dir> <sample_info> <in_seurat_obj> <in_matrix_dir> <out_dir> <prefix>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
SEURAT_SCRIPT_DIR=$1
SAMPLE_INFO=$2
IN_SEURAT_OBJ=$3
IN_MATRIX_DIR=$4
OUT_DIR=$5
PREFIX=$6

OBS="Sample name,Diagnosis,Region code"

SEURAT_DIR="${OUT_DIR}/seurat"
ADATA_DIR="${OUT_DIR}/adata"
MATRIX_DIR="${OUT_DIR}/matrix/${PREFIX}"

[ -d "${SEURAT_DIR}" ] || mkdir -p "${SEURAT_DIR}"
[ -d "${ADATA_DIR}" ] || mkdir -p "${ADATA_DIR}"
[ -d "${MATRIX_DIR}" ] || mkdir -p "${MATRIX_DIR}"

OUT_SEURAT_OBJ="${SEURAT_DIR}/${PREFIX}.Rds"
META_FILE="${SEURAT_DIR}/${PREFIX}_metadata.tsv"
OUT_ADATA_OBJ="${ADATA_DIR}/${PREFIX}.h5ad"

FEATURES="${IN_MATRIX_DIR}/features.tsv.gz"
FEATURE_MAPPING="${MATRIX_DIR}/feature_mapping.tsv"

Rscript ${SCRIPT_DIR}/rename_seurat_meta_for_GCA.R "${IN_SEURAT_OBJ}" "${OUT_SEURAT_OBJ}" "${SAMPLE_INFO}"
Rscript ${SEURAT_SCRIPT_DIR}/generate_metadata.R ${OUT_SEURAT_OBJ} ${META_FILE}
Rscript ${SEURAT_SCRIPT_DIR}/generate_matrix.R ${OUT_SEURAT_OBJ} ${MATRIX_DIR}
Rscript ${SEURAT_SCRIPT_DIR}/feature_mapping.R "${FEATURES}" "${FEATURE_MAPPING}"
python ${SCRIPT_DIR}/create_scampy_object.py "${MATRIX_DIR}" "${META_FILE}" "${FEATURE_MAPPING}" "${OUT_ADATA_OBJ}"
python ${SCRIPT_DIR}/subset_obs.py "${OUT_ADATA_OBJ}" "${OBS}" "${OUT_ADATA_OBJ}"

exit

