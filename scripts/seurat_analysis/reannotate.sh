
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt 6 ]]
then
    echo "Usage: bash $0 <obj> <cl_mode> <file_list> <annotation> <plot_prefix> <reduction>"
fi

OBJECT=$1
CL_MODE=$2
FILE_LIST=$3
ANNOTATION=$4
PLOT_PREFIX=$5
REDUCTION=$6

OUT_CL_MODE="${CL_MODE}_subclusters"

Rscript ${SCRIPT_DIR}/add_subclusters.R "${OBJECT}" ${CL_MODE} "${FILE_LIST}" ${OUT_CL_MODE}
Rscript ${SCRIPT_DIR}/cluster_labels.R "${OBJECT}" ${OUT_CL_MODE} "${ANNOTATION}" "${PLOT_PREFIX}" ${REDUCTION}

exit

