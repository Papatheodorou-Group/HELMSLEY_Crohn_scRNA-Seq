PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
R_SCRIPT_DIR="../../../r_scripts"
GIT_DIR="${R_SCRIPT_DIR}/gene_expression_analysis"

EXP_N=( "004_N_E_S" \
        "006_N_E_S" \
        "007_N_E_S" \
        "012_N_E_S" )
EXP_C=( "003_C_E_S" \
        "008_C_E_S" \
        "010_C_E_S" \
        "013_C_E_S" )
EXP=( ${EXP_N[@]} ${EXP_C[@]} ) 

CL_MODE="clusters_pca_vst_top5000_k30"
RES=( $(cat "param_clust.txt" | grep '^res=' | sed "s/res=//g" | tr ',' ' ') )

MEM="30"

DEA="4_DEA.sh"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

batch_job () {

    EXP_ID=$1
    SET=$2
    SLOT=$3
    SCRIPT=$4

    OUT_OBJ="${EXP_ID}/${SET}"
    
    # Run DEA 
    ID="DEA_${EXP_ID}"
    COMMAND="${SCRIPT_DIR}/${SCRIPT} \"${GIT_DIR}\" \"${OUT_OBJ}\" ${SLOT}" 
    sbatch -t 24:00:00 --mem=${MEM}G \
           -J ${ID} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"
}

# run differential expression analysis

i=0
for exp in ${EXP[@]}
do  
    for res in ${RES[@]}
    do 
        SLOT="${CL_MODE}_res${res}"
        batch_job ${exp} ${exp}_ccRegress ${SLOT} ${DEA}
    done
done


exit

