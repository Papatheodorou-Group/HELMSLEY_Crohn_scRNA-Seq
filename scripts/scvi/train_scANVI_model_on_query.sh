source ~/.bashrc
conda activate scvi-env

# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

META=( "category" "Integrated_05" ) # broad or fine-grained cell-type annotation

if [[ $# -lt "4" ]]
then
    echo "usage: bash create_model_with_scANVI.sh <q_obj_prefix> <r_mod_out> <reduct_name_an> <stat_prefix> <GCA_stats>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
Q_OBJ_PREFIX=$1
R_MOD_OUT_PREFIX=$2
REDUCT_NAME_AN=$3
STAT_PREFIX=$4
GCA_STATS=$5

for m in ${META[@]}
do
    R_MOD_OUT="${R_MOD_OUT_PREFIX}_$m" 
    python ${SCRIPT_DIR}/create_scANVI_model_query.py "${Q_OBJ_PREFIX}" "${R_MOD_OUT}" ${m} ${REDUCT_NAME_AN}
    # compute soft label assignment statistics
    Rscript ${SCRIPT_DIR}/fuzzy_cell_type_matrix.R "${Q_OBJ_PREFIX}_soft${m}.tsv" "${GCA_STATS}" "${STAT_PREFIX}_soft${m}_fuzzyMatrix"
done

Rscript ${SCRIPT_DIR}/pc_soft_annotation.R "${Q_OBJ_PREFIX}_soft${META[1]}.tsv" "${GCA_STATS}" "${Q_OBJ_PREFIX}_soft${META[1]}"
Rscript ${SCRIPT_DIR}/plot_fuzzy_cell_type_matrix.R "${STAT_PREFIX}_soft${META[1]}_fuzzyMatrix"

# query
python ${SCRIPT_DIR}/print_obs.py "${Q_OBJ_PREFIX}.h5ad" "${STAT_PREFIX}.tsv"
Rscript ${SCRIPT_DIR}/stats_barplot.R "${STAT_PREFIX}.tsv" "${STAT_PREFIX}"

exit

