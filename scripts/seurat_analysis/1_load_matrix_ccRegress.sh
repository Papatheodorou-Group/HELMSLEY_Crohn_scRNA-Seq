
MIN_GENES="0"
MIN_CELLS="3"

if [[ $# -lt "3" ]]
then
    echo "bash $0 <script_dir> <matrix_list> <in_dir> <out_dir>"
    exit
fi

SCRIPT_DIR="$1"
MATRIX_LIST="$2"
IN_DIR="$3"
OUT_DIR="$4"

R_SCRIPT_LOAD_MATRIX="${SCRIPT_DIR}/Seurat_load_multi_matrix.R"

# create the object
echo "create the Seurat object"
Rscript "${R_SCRIPT_LOAD_MATRIX}" \
        --matrix_dir_list "${MATRIX_LIST}" \
        --out_RDS "${OUT_DIR}/object.Rds" \
        --min_genes ${MIN_GENES} \
        --min_cells ${MIN_CELLS} \
        --vars_to_regress "S.Score,G2M.Score" \
        1> "${OUT_DIR}/Seurat_load_multi_matrix.STDOUT" \
        2> "${OUT_DIR}/Seurat_load_multi_matrix.STDERR"



exit


