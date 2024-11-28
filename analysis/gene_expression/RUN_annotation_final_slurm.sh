
## N.B.!!!
## A manual annotation file is built after RUN_subclustering.sh,
## called ANNOTATION_merged_N_filt_2ndRound_TIL.tsv for normal (with doublet removal before annotation)
## and ANNOTATION_merged_C_filt_no010.tsv for disease (without doublet removal; cluster filtering at the annotation stage)

#### PARAMETERS #####

PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
WDIR=$(pwd)

DR_MODE="pca_vst_top5000"
CL_MODE="clusters_${DR_MODE}_k30_res1.2"

MEM="50"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}



### EXECUTE ###


## normal

DIR="merged/merged_N_filt_2ndRound_TIL"

PLOT_PREFIX="${DIR}/plots/${DR_MODE}/UMAP_plot_${CL_MODE}"

OBJECT="${DIR}/hvg_pca_clust/object.Rds"
FILE_LIST="${DIR}/FILE_LIST.tsv"
ANNOTATION="ANNOTATION_merged_N_filt_2ndRound_TIL.tsv"

ID="reannotate_N"
COMMAND="${SCRIPT_DIR}/reannotate.sh \"${OBJECT}\" ${CL_MODE} \"${FILE_LIST}\" \"${ANNOTATION}\" \"${PLOT_PREFIX}\" \"umap_${DR_MODE}\"" 
sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"


## crohn

DIR="merged/merged_C_filt_no010"

PLOT_PREFIX="${DIR}/plots/${DR_MODE}/UMAP_plot_${CL_MODE}"

OBJECT="${DIR}/hvg_pca_clust/object.Rds"
FILE_LIST="${DIR}/FILE_LIST.tsv"
ANNOTATION="ANNOTATION_merged_C_filt_no010.tsv"

ID="reannotate_C"
COMMAND="${SCRIPT_DIR}/reannotate.sh \"${OBJECT}\" ${CL_MODE} \"${FILE_LIST}\" \"${ANNOTATION}\" \"${PLOT_PREFIX}\" \"umap_${DR_MODE}\"" 
sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"

exit

