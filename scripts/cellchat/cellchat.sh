
SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt "3" ]]
then
    echo ""
    echo "Usage: bash $0 <OBJ1> <OBJ2> <OUT_DIR>"
    echo ""
    exit
fi

OBJ1="$1"
OBJ2="$2"
OUT_DIR="$3"

Rscript "${SCRIPT_DIR}/cellchat_compute.R" "${OBJ1}" "${OBJ2}" "${OUT_DIR}"
Rscript "${SCRIPT_DIR}/cellchat_plot.R" "${OUT_DIR}"


exit


