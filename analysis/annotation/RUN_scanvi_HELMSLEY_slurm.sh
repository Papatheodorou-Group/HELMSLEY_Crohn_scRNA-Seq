PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR="../../scripts/scvi"
SEURAT_SCRIPT_DIR="../../scripts/seurat_analysis"

SEURAT_DIR="../gene_expression/merged"
MATRIX_DIR="../cell_calling"
Q_DATA_DIR="../../data"

MEM="30"
CORES="2"

############### GCA ################

OUT_AN_MODEL_DIR="GCA/scanvi_models"
OUT_OBJ_DIR="GCA/adata"
OUT_STAT_DIR="GCA/stats"
GCA_STATS="GCA/GCA_stats.tsv"

SET="adult_pediatric"

SCANVI_REDUCT="X_scANVI_category" # if weight_decay = 0.0 works as expected, the reference coordinates should not change upon training the query

R_OBJ="${OUT_OBJ_DIR}/${SET}.h5ad"
  
############ HELMSLEY ##############
       
OUT_Q_DIR="HELMSLEY"

# Q_SETS=( "merged_N_filt_TIL" \
#         "merged_N_filt_2ndRound_TIL" \
#         "merged_C_filt" \
#         "merged_C_filt_2ndRound" \
#         "merged_C_filt_no010" \
#         "merged_C_filt_2ndRound_no010" \
#         "merged_all_filt" \
#         "merged_all_filt_2ndRound" )
Q_SETS=( "merged_all_filt" \
         "merged_all_filt_2ndRound" )
IN_QUERY_OBJECTS=( $( for exp in ${Q_SETS[@]} ; do echo "${SEURAT_DIR}/${exp}/object.Rds" ; done ) )

IN_QUERY_MATRIX="${MATRIX_DIR}/SR1803_20220727_002_N_S_S/outs/filtered_feature_bc_matrix" # this is only needed for ENSEMBL gene IDs

SAMPLE_INFO="${Q_DATA_DIR}/scRNA_samples_data.tsv"

Q_PLOT_DIR="${OUT_Q_DIR}/plots"
Q_STAT_DIR="${OUT_Q_DIR}/stats"
Q_STAT_DIR_FILT="${OUT_Q_DIR}/stats_filt"

[ -d "${OUT_STAT_DIR}" ] || mkdir -p "${OUT_STAT_DIR}"
STAT_R_PREFIX="${OUT_STAT_DIR}/${SET}"

#####################################

# map the HELMSLEY data to the GCA model

SEURAT2ADATA="${SCRIPT_DIR}/seurat_to_adata.sh"
TRAIN="${SCRIPT_DIR}/train_scANVI_model_on_query.sh"
PLOT="${SCRIPT_DIR}/plot_celltypes_prediction.sh"
FILTER="${SCRIPT_DIR}/filter_cells_by_annotation_consistency.sh"

MODEL_AN="${OUT_AN_MODEL_DIR}/${SET}"
STAT_R="${STAT_R_PREFIX}.tsv" # this is just to retrieve the cell types from GCA

[ -d "${Q_STAT_DIR}" ] || mkdir -p "${Q_STAT_DIR}"
[ -d "${Q_STAT_DIR_FILT}" ] || mkdir -p "${Q_STAT_DIR_FILT}"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

M=${#Q_SETS[@]}
for (( j=0; j<$M; j++ ))
do
    # convert the query object from Seurat to anndata
    ID3="H_${j}_step3"
    COMMAND="${SEURAT2ADATA} \"${SEURAT_SCRIPT_DIR}\" \"${SAMPLE_INFO}\" \"${IN_QUERY_OBJECTS[$j]}\" \
                             \"${IN_QUERY_MATRIX}\" \"${OUT_Q_DIR}\" \"${Q_SETS[$j]}\""
    J=$(sbatch -t 10:00:00 --mem=${MEM}G \
               -J ${ID3} \
               -o ${OE}/${ID3}.STDOUT \
               -e ${OE}/${ID3}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")
    J3=${J#Submitted batch job }

    Q_OBJ_PREFIX="${OUT_Q_DIR}/adata/${Q_SETS[$j]}"
    Q_OBJ="${Q_OBJ_PREFIX}.h5ad"
    STAT_Q_PREFIX="${Q_STAT_DIR}/${Q_SETS[$j]}"
 
    # map the query cells to the reference model (2 models used, one per "SET")
    # perform both soft and hard annotation
    ID4="H_${j}_step4"
    COMMAND="${TRAIN} \"${Q_OBJ_PREFIX}\" \"${MODEL_AN}\" scANVI \"${STAT_Q_PREFIX}\" \"${GCA_STATS}\""
    J=$(sbatch -t 10:00:00 --mem=${MEM}G -c ${CORES} \
               -J ${ID4} \
               --dependency=afterok:${J3} \
               -o ${OE}/${ID4}.STDOUT \
               -e ${OE}/${ID4}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")
    J4=${J#Submitted batch job }

    # directory containing the plot for a specific query set
    OUT_Q_PLOT_DIR="${Q_PLOT_DIR}/${Q_SETS[$j]}"
    [ -d "${OUT_Q_PLOT_DIR}" ] || mkdir -p "${OUT_Q_PLOT_DIR}"

    # plot cell types (all)
    ID5="H_${j}_step5"
    COMMAND="${PLOT} \"${Q_OBJ}\" \"${R_OBJ}\" ${SCANVI_REDUCT} \"${OUT_Q_PLOT_DIR}\""
    sbatch -t 10:00:00 --mem=${MEM}G -c ${CORES} \
           -J ${ID5} \
           --dependency=afterok:${J4} \
           -o ${OE}/${ID5}.STDOUT \
           -e ${OE}/${ID5}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"
    
    STAT_Q="${STAT_Q_PREFIX}.tsv"
    STAT_Q_PREFIX_FILT="${Q_STAT_DIR_FILT}/${Q_SETS[$j]}"
    Q_OBJ_FILT="${OUT_Q_DIR}/adata/${Q_SETS[$j]}_filt.h5ad"
    
    # select consistent cells: the category deduced from cell-type prediction and the predicted category are the same
    ID6="H_${j}_step6"
    COMMAND="${FILTER} \"${STAT_Q}\" \"${STAT_R}\" \"${Q_OBJ}\" \"${Q_OBJ_FILT}\" \"${STAT_Q_PREFIX_FILT}\""
    J=$(sbatch -t 10:00:00 --mem=${MEM}G \
               -J ${ID6} \
               --dependency=afterok:${J4} \
               -o ${OE}/${ID6}.STDOUT \
               -e ${OE}/${ID6}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")
    J6=${J#Submitted batch job }

    # directory containing the plot for a specific query set
    OUT_Q_PLOT_DIR_FILT="${Q_PLOT_DIR}/${Q_SETS[$j]}_filtByAnn"
    [ -d "${OUT_Q_PLOT_DIR_FILT}" ] || mkdir -p "${OUT_Q_PLOT_DIR_FILT}"

    # plot cell types (filt)
    ID7="H_${j}_step7"
    COMMAND="${PLOT} \"${Q_OBJ_FILT}\" \"${R_OBJ}\" ${SCANVI_REDUCT} \"${OUT_Q_PLOT_DIR_FILT}\""
    sbatch -t 10:00:00 --mem=${MEM}G \
           -J ${ID7} \
           --dependency=afterok:${J6} \
           -o ${OE}/${ID7}.STDOUT \
           -e ${OE}/${ID7}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"

done

exit


