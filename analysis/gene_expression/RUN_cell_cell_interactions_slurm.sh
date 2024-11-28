
## cell-cell interaction analysis depends on annotation
## there are 2 levels od annotation: cell_type and category
## I am using "category" here, as reported in the 3rd column of "ANNOTATION_merged_*" 
## This info is stored in "category_final" in the Seurat object
## here, remove 'NA' cells (label convention for unclassified cells) and labels that 
## are present in only one condition (here, "Epithelial inflamed" from Crohn's).

#### PARAMETERS #####

PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/cellchat" ; pwd )
WDIR=$(pwd)

MEM="50"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}



### EXECUTE ###

OBJECT_N="merged/merged_N_filt_2ndRound_TIL/hvg_pca_clust/object.Rds"
OBJECT_C="merged/merged_C_filt_no010/hvg_pca_clust/object.Rds"
OUT_DIR="cellchat"

ID="cellchat"
COMMAND="${SCRIPT_DIR}/cellchat.sh \"${OBJECT_N}\" \"${OBJECT_C}\" \"${OUT_DIR}\"" 
sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"



exit

