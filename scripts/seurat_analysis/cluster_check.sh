
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt 5 ]]
then
    echo "Usage: bash $0 <obj> <cl_mode> <dea> <tsv> <plot>"
fi

OBJECT=$1
CL_MODE=$2
DEA=$3
TSV=$4
PLOT=$5

Rscript ${SCRIPT_DIR}/cluster_check.R "${OBJECT}" ${CL_MODE} "${DEA}" "${TSV}" "${PLOT}"

exit

