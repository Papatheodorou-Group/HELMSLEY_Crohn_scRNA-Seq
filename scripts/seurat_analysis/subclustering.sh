
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt 7 ]]
then
    echo "Usage: bash $0 <obj> <cl_mode> <cl_id> <res> <cat> <plot_prefix> <dea>"
    exit
fi

OBJECT=$1
CL_MODE=$2
CL_ID=$3
RES=$4
CAT=$5
PLOT_PREFIX=$6
DEA=$7

Rscript ${SCRIPT_DIR}/subclustering.R "${OBJECT}" ${CL_MODE} ${CL_ID} "${RES}" "${CAT}" "${PLOT_PREFIX}" "${DEA}"

exit

