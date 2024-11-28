
#### PARAMETERS #####

PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR="../../../r_scripts/gene_expression_analysis"
WDIR=$(pwd)

DR_MODE="pca_vst_top5000"
CL_MODE="clusters_${DR_MODE}_k30_res1.2"

MEM="60"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}



#### FUNCTIONS ####

# to decide whether a cluster should be removed or not
# (high expression of cell death genes? co-expression of incompatible markers?)  
cl_check () {

    LABEL=$1
    DEA=$2
    TSV=$3
    PLOT=$4
        
    # extract the top markers and plot their average expression across cell types
    ID="cluster_check_${LABEL}"
    COMMAND="${SCRIPT_DIR}/cluster_check.sh \"${OBJECT}\" ${CL_MODE} \"${DEA}\" \"${TSV}\" \"${PLOT}\""
    J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
}

# recluster the cells 
# N.B. do not need dependency anymore, as the object won't be updated -> save metadata to tsv file instead
subclust () {

    LABEL=$1
    OBJECT=$2
    CL_SLOT=$3
    CL=$4
    RES=$5
    CAT=$6
    PLOT=$7
    DEA=$8
    
    # recluster a specific cluster
    ID="subcluster_${LABEL}"
    COMMAND="${SCRIPT_DIR}/subclustering.sh \"${OBJECT}\" ${CL_SLOT} ${CL} ${RES} \"${CAT}\" \"${PLOT}\" \"${DEA}\"" 
    J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")

    echo $J
}



### EXECUTE ###


## normal

DIR="merged/merged_N_filt_2ndRound_TIL"

DEA_DIR="${DIR}/DEA/${CL_MODE}/MAST"
PLOT_PREFIX="${DIR}/plots/${DR_MODE}/UMAP_plot_${CL_MODE}"

OBJECT="${DIR}/hvg_pca_clust/object.Rds"
FILE_LIST="${DIR}/FILE_LIST.tsv"

CL_CHECK=( "1" "2" "4" "6" "7" "11" "12" "17" "18" "19" "21" "28" )

SUBCLUST=( "3" "13" "22" "14" "15" "23" "20" "26" )
RES=( "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" )
CAT=( 'T cells' 'T cells' 'T cells' 'B cells' 'B cells' 'B cells' 'Myeloid' 'Mesenchymal' )

for cl in ${CL_CHECK[@]}
do
    dea="${DEA_DIR}/DEG_MAST_cl${cl}-all.tsv"
    tsv="${DEA_DIR}/DEG_MAST_cl${cl}_AvgExpTopDEGs.tsv"
    plot="${DEA_DIR}/DEG_MAST_cl${cl}_AvgExpTopDEGs.pdf"

#    cl_check ${cl}_normal ${dea} ${tsv} ${plot}
done

cp /dev/null ${FILE_LIST}
dep=""
for (( i=0; i<${#SUBCLUST[@]}; i++ ))
do
    plot="${PLOT_PREFIX}_cl${SUBCLUST[$i]}_subcl"
    dea="${DEA_DIR}/DEG_MAST_cl${SUBCLUST[$i]}_subclAllMarkers.tsv"

#    subclust ${SUBCLUST[$i]}_normal ${OBJECT} ${CL_MODE} ${SUBCLUST[$i]} ${RES[$i]} "${CAT[$i]}" "${plot}" "${dea}" 
    echo "${SUBCLUST[$i]}	${plot}_clusters.tsv" >> ${FILE_LIST} 
done


## crohn

DIR="merged/merged_C_filt_no010"

DEA_DIR="${DIR}/DEA/${CL_MODE}/MAST"
PLOT_PREFIX="${DIR}/plots/${DR_MODE}/UMAP_plot_${CL_MODE}"

OBJECT="${DIR}/hvg_pca_clust/object.Rds"
FILE_LIST="${DIR}/FILE_LIST.tsv"

CL_CHECK=( "1" "5" "7" "9" "11" "17" "18" "24" )

SUBCLUST=( "2" "3" "4" "10" "14" "19" "20" "21" "23" "29" )
RES=( "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" "0.3" )
CAT=( 'T cells' 'T cells' 'B cells' 'Myeloid' 'Mesenchymal' 'Mesenchymal' 'Mesenchymal' 'B cells' 'Mesenchymal' )

for cl in ${CL_CHECK[@]}
do
    dea="${DEA_DIR}/DEG_MAST_cl${cl}-all.tsv"
    tsv="${DEA_DIR}/DEG_MAST_cl${cl}_AvgExpTopDEGs.tsv"
    plot="${DEA_DIR}/DEG_MAST_cl${cl}_AvgExpTopDEGs.pdf"

#    cl_check ${cl}_crohn ${dea} ${tsv} ${plot}
done

cp /dev/null ${FILE_LIST}
dep=""
for (( i=0; i<${#SUBCLUST[@]}; i++ ))
do
    plot="${PLOT_PREFIX}_cl${SUBCLUST[$i]}_subcl"
    dea="${DEA_DIR}/DEG_MAST_cl${SUBCLUST[$i]}_subclAllMarkers.tsv"

#    subclust ${SUBCLUST[$i]}_crohn ${OBJECT} ${CL_MODE} ${SUBCLUST[$i]} ${RES[$i]} "${CAT[$i]}" "${plot}" "${dea}"

    echo "${SUBCLUST[$i]}	${plot}_clusters.tsv" >> ${FILE_LIST}
done


exit

