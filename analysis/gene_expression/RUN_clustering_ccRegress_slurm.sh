PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR="../../../r_scripts/gene_expression_analysis" # this path is for the new R setup (multiprocess -> multicore)
WDIR=$(pwd)

EXP=( "002_N_S_S" \
      "003_C_E_S" \
      "004_N_E_S" \
      "006_N_E_S" \
      "007_N_E_S" \
      "008_C_E_S" \
      "010_C_E_S" \
      "012_N_E_S" \
      "013_C_E_S" )

MEM="25"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

batch_job () {

    SET=$1
    LABEL=$2

    FILE_LIST="matrix_dir_${LABEL}.tsv"

    OUT_OBJ="${LABEL}/${SET}_ccRegress"
    [[ -d "${OUT_OBJ}" ]] || mkdir -p "${OUT_OBJ}"

    # Generate the seurat object 
    COMMAND="${SCRIPT_DIR}/1_load_matrix_ccRegress.sh \"${GIT_DIR}\" \"${FILE_LIST}\" \"${LABEL}\" \"${OUT_OBJ}\"" 
    J=$(sbatch -t 24:00:00 --mem=${MEM}G \
               -J step1_${LABEL} \
               -o ${OE}/step1_${LABEL}.STDOUT \
               -e ${OE}/step1_${LABEL}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")
    J1=${J#Submitted batch job }

    # Run the clustering
    # error with ggplot in plot script!!! 
    # FIXED: remotes::install_github("thomasp85/patchwork")
    #        remotes::install_github("tidyverse/ggplot2", ref = remotes::github_pull("5592"))
    COMMAND="${SCRIPT_DIR}/2_clustering.sh \"${GIT_DIR}\" \"${OUT_OBJ}\" \"${WDIR}\"" 
    J=$(sbatch -t 24:00:00 --mem=${MEM}G \
               -J step2_${LABEL} \
               --dependency=afterok:${J1} \
               -o ${OE}/step2_${LABEL}.STDOUT \
               -e ${OE}/step2_${LABEL}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")
    J2=${J#Submitted batch job }
       
    # collect statistics on generated clustering solutions
    COMMAND="${SCRIPT_DIR}/3_collect_cl_stat.sh \"${WDIR}\" \"${OUT_OBJ}\"" 
    sbatch -t 24:00:00 --mem=${MEM}G \
           -J step3_${LABEL} \
           --dependency=afterok:${J2} \
           -o ${OE}/step3_${LABEL}.STDOUT \
           -e ${OE}/step3_${LABEL}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"
}

for exp in ${EXP[@]}
do
    FILE_LIST="matrix_dir_${exp}.tsv"
    batch_job ${exp} ${exp}
    
done


exit

