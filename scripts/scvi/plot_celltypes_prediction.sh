source ~/.bashrc
conda activate scvi-env

META=( "category" "Integrated_05" ) # broad or fine-grained cell-type annotation

if [[ $# -lt "4" ]]
then
    echo "usage: bash plot_cell_type_prediction.sh <q_obj> <r_obj> <reduct_name> <plot_dir>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
Q_OBJ=$1
R_OBJ=$2
REDUCT_NAME=$3
PLOT_DIR=$4

python ${SCRIPT_DIR}/plot_celltypes_prediction.py "${Q_OBJ}" "${R_OBJ}" "${PLOT_DIR}" ${REDUCT_NAME} ${META[@]}

# MISC_SCRIPT_DIR=$( cd "${SCRIPT_DIR}/../misc"; pwd )
# FIXME: missing poppler on codon
# cd "${PLOT_DIR}"
# bash ${MISC_SCRIPT_DIR}/pdf2png.sh

exit

