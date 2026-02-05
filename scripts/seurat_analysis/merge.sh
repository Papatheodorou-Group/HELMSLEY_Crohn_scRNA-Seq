
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )
# DATA_DIR=$( cd ${SCRIPT_DIR}; cd "../../data" ; pwd )

# GENES="${DATA_DIR}/Cells_genes_min.tsv"

if [[ $# -lt 3 ]]
then
    echo "Usage: bash $0 <dir> <dir1> ... <dirn>"
    exit
fi

DIR=$1
IN_DIR=${@:2:$#}
# SUFFIX=$2
# IN_DIR=${@:3:$#}

[[ -d "${DIR}" ]] || mkdir -p "${DIR}"

IN_OBJ=( )
for in_dir in ${IN_DIR[@]}
do
#    IN_OBJ=( ${IN_OBJ[@]} "${in_dir}/${in_dir}_${SUFFIX}/object.Rds" )
    IN_OBJ=( ${IN_OBJ[@]} "${in_dir}/object.Rds" )
done

# merge objects
Rscript ${SCRIPT_DIR}/merge_obj.R "${DIR}" ${IN_OBJ[@]}

# [ -d "${DIR}/umaps_genes" ] || mkdir -p "${DIR}/umaps_genes"

# plot signature genes
# Rscript ${SCRIPT_DIR}/plot_genes.R "${DIR}/umaps_genes" "${DIR}/object.Rds" ${GENES}

exit

