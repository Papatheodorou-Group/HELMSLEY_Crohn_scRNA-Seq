
PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR="../../scripts/scvi"
SEURAT_SCRIPT_DIR="../../scripts/seurat_analysis"
DATA_DIR="../../data/GCA"

MEM="50"
CORES="2"

############### GCA ################

OUT_MODEL_DIR="GCA/scvi_models"
OUT_AN_MODEL_DIR="GCA/scanvi_models"
OUT_OBJ_DIR="GCA/adata"
OUT_PLOT_DIR="GCA/plots"
OUT_STAT_DIR="GCA/stats"

IN_OBJ="${DATA_DIR}/obj_healthy_adult_pediatric.h5ad"

SET="adult_pediatric"

SCVI_REDUCT="scVI"
SCANVI_REDUCT="scANVI"

META=( "category" "Integrated_05" )

#####################################

[ -d "${OUT_MODEL_DIR}" ] || mkdir -p "${OUT_MODEL_DIR}"
[ -d "${OUT_AN_MODEL_DIR}" ] || mkdir -p "${OUT_AN_MODEL_DIR}"
[ -d "${OUT_OBJ_DIR}" ] || mkdir -p "${OUT_OBJ_DIR}"
[ -d "${OUT_STAT_DIR}" ] || mkdir -p "${OUT_STAT_DIR}"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

# generate scVI and scANVI models on GCA

MODEL="${OUT_MODEL_DIR}/${SET}"
MODEL_AN="${OUT_AN_MODEL_DIR}/${SET}"
OUT_OBJ="${OUT_OBJ_DIR}/${SET}"
STAT_R_PREFIX="${OUT_STAT_DIR}/${SET}"
OUT_STATS="GCA/GCA_stats.tsv"

OUT_DIR="${OUT_PLOT_DIR}/${SET}"
[ -d "${OUT_DIR}" ] || mkdir -p "${OUT_DIR}"

R_OBJ="${OUT_OBJ}.h5ad"

# train the models and save them on anndata structure
ID1="step1"
COMMAND="${SCRIPT_DIR}/create_models_on_ref.sh \"${IN_OBJ}\" \"${MODEL}\" \"${MODEL_AN}\" \"${OUT_OBJ}\" ${SCVI_REDUCT} ${SCANVI_REDUCT} \"${STAT_R_PREFIX}\""
J=$(sbatch -t 48:00:00 --mem=${MEM}G -c ${CORES} \
           -J ${ID1} \
           -o ${OE}/${ID1}.STDOUT \
           -e ${OE}/${ID1}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J1=${J#Submitted batch job }

# plot the cells on PCA or scvi latent space
ID2="step2"
COMMAND="${SCRIPT_DIR}/plot_celltypes.sh \"${R_OBJ}\" \"${OUT_DIR}\" X_${SCVI_REDUCT}"
J=$(sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID2} \
           --dependency=afterok:${J1} \
           -o ${OE}/${ID2}.STDOUT \
           -e ${OE}/${ID2}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J2=${J#Submitted batch job }  

# plot the cells on PCA or scanvi latent space
for meta in ${META[@]}
do
    ID3="step3_${meta}"
    COMMAND="${SCRIPT_DIR}/plot_celltypes.sh \"${R_OBJ}\" \"${OUT_DIR}\" X_${SCANVI_REDUCT}_${meta}"
    sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID3} \
           --dependency=afterok:${J1} \
           -o ${OE}/${ID3}.STDOUT \
           -e ${OE}/${ID3}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"
done

# collect stats
ID4="step4"
COMMAND="${SCRIPT_DIR}/collect_stats.sh \"${OUT_OBJ}\" \"${OUT_STATS}\""
sbatch -t 2:00:00 --mem=${MEM}G \
       -J ${ID4} \
       --dependency=afterok:${J1} \
       -o ${OE}/${ID4}.STDOUT \
       -e ${OE}/${ID4}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"

exit

