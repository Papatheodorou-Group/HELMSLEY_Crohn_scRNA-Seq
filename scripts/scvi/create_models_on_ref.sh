source ~/.bashrc
conda activate scvi-env

# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

META=( "category" "Integrated_05" ) # broad or fine-grained cell-type annotation

if [[ $# -lt "7" ]]
then
    echo "usage: bash create_scVI_model.sh <in_obj> <mod> <mod_an> <out_obj> <reduct_name> <reduct_name_an> <stat_prefix>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
IN_OBJ=$1
MOD=$2
MOD_AN=$3
OUT_OBJ=$4
REDUCT_NAME=$5
REDUCT_NAME_AN=$6
STAT_PREFIX=$7

OBS="Sample name,Diagnosis,Region code,category,Integrated_05"

python ${SCRIPT_DIR}/subset_obs.py "${IN_OBJ}" "${OBS}" "${OUT_OBJ}.h5ad"

# reference
python ${SCRIPT_DIR}/print_obs.py "${OUT_OBJ}.h5ad" "${STAT_PREFIX}.tsv"
Rscript ${SCRIPT_DIR}/stats_barplot.R "${STAT_PREFIX}.tsv" "${STAT_PREFIX}"

python ${SCRIPT_DIR}/create_scVI_model.py "${OUT_OBJ}.h5ad" "${MOD}" "${OUT_OBJ}.h5ad" ${REDUCT_NAME}

for m in ${META[@]}
do
    R_MOD_OUT="${MOD_AN}_$m" 
    python ${SCRIPT_DIR}/create_scANVI_model_ref.py "${OUT_OBJ}" "${MOD}" "${R_MOD_OUT}" ${m} ${REDUCT_NAME_AN}
done 

exit

