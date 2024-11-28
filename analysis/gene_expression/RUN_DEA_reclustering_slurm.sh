
PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR="../../../r_scripts/gene_expression_analysis"

EXP="merged"
EXP_TYPE=( "merged_N_filt_TIL" \
           "merged_C_filt_no010" \
           "merged_N_filt_2ndRound_TIL" \
           "merged_C_filt_2ndRound_no010" )

CL_MODE="clusters_pca_vst_top5000_k30"
RES=( $(cat "param_clust_filt.txt" | grep '^res=' | sed "s/res=//g" | tr ',' ' ') )

MEM="50"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

DEA="4_DEA.sh"

batch_job () {

    EXP_ID=$1
    SET=$2
    SLOT=$3
    SCRIPT=$4

    OUT_OBJ="${EXP_ID}/${SET}"
    
    # Run DEA 
    ID="DEA_${SET}"
    COMMAND="${SCRIPT_DIR}/${SCRIPT} \"${GIT_DIR}\" \"${OUT_OBJ}\" ${SLOT}" 
    sbatch -t 30:00:00 --mem=${MEM}G \
              -J ${ID} \
              -o ${OE}/${ID}.STDOUT \
              -e ${OE}/${ID}.STDERR \
              --wrap="${PREAMBLE} ; bash ${COMMAND}"

}

# run differential expression analysis
# for (( i=0; i < ${#CL_MODE[@]}; i++ ))
# do  
#    batch_job ${EXP} ${EXP_TYPE[$i]} ${CL_MODE[$i]} ${DEA}
# done

# NEW: test all resolutions
for exp in ${EXP_TYPE[@]}
do  
    for res in ${RES[@]}
    do 
        SLOT="${CL_MODE}_res${res}"
        batch_job ${EXP} ${exp} ${SLOT} ${DEA}
    done
done


exit

