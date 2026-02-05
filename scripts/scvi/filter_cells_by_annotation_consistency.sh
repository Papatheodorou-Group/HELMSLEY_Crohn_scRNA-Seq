source ~/.bashrc
conda activate scvi-env

# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

if [[ $# -lt "5" ]]
then
    echo "usage: bash filter_cells_by_annotation_consistency.sh <stat_q> <stat_r> <in_q_obj> <out_q_obj> <stat_prefix>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
STAT_Q=$1
STAT_R=$2
IN_Q_OBJ=$3
OUT_Q_OBJ=$4
STAT_PREFIX=$5

CELLS="$( dirname ${STAT_Q} )/selected_cells.txt"  

Rscript ${SCRIPT_DIR}/select_consistent_cells.R "${STAT_Q}" "${STAT_R}" "${CELLS}"
python ${SCRIPT_DIR}/filter_obj_cells.py "${IN_Q_OBJ}" "${CELLS}" "${OUT_Q_OBJ}"

# query
python ${SCRIPT_DIR}/print_obs.py "${OUT_Q_OBJ}" "${STAT_PREFIX}.tsv"
Rscript ${SCRIPT_DIR}/stats_barplot.R "${STAT_PREFIX}.tsv" "${STAT_PREFIX}"

exit

