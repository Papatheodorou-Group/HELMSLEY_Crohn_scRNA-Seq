PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
GIT_DIR="../../../r_scripts/gene_expression_analysis"
WDIR=$(pwd)

EXP_N=( "004_N_E_S" \
        "006_N_E_S" \
        "007_N_E_S" \
        "012_N_E_S" )
EXP_C=( "003_C_E_S" \
        "008_C_E_S" \
        "010_C_E_S" \
        "013_C_E_S" )
EXP_C_no010=( "003_C_E_S" \
              "008_C_E_S" \
              "013_C_E_S" )
EXP=( ${EXP_N[@]} ${EXP_C[@]} )      

SUFFIX="ccRegress_filt_2ndRound"
IN_DIR_N=( $( for exp in ${EXP_N[@]} ; do echo ${exp}/${exp}_${SUFFIX} ; done ) )
IN_DIR_C=( $( for exp in ${EXP_C[@]} ; do echo ${exp}/${exp}_${SUFFIX} ; done ) )
IN_DIR_C_no010=( $( for exp in ${EXP_C_no010[@]} ; do echo ${exp}/${exp}_${SUFFIX} ; done ) )
IN_DIR=( ${IN_DIR_N[@]} ${IN_DIR_C_no010[@]} )

# select the clustering solution for 2nd filtering round
CL_MODE=( "clusters_pca_vst_top5000_k30_res0.6" \
          "" \
          "" \
          "clusters_pca_vst_top5000_k30_res1.2" 
          "" \
          "clusters_pca_vst_top5000_k30_res0.8" \
          "clusters_pca_vst_top5000_k30_res0.7" \
          "" ) 
          
# cluster ID -> to be removed, as they contain circa > 70% doublets       
CL_FILT=( "13" "" "" "21" "" "22,23" "20" "" ) 

# dimensional reduction where doublets have been calculated
DR_DOUBL="pca_vst_top5000" 
    
MEM="20"
MEM2="30"

FILT="cell_filtering_2ndRound.sh"
MERGE="merge.sh"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

batch_job () {

    EXP_ID=$1
    SET=$2
    DR=$3
    SLOT=$4
    CL=$5

    DIR_CLUSTERS="${EXP_ID}/${SET}" # dir where clustering is stored
    DIR_DOUBLETS="${EXP_ID}/${SET}_filt" # dir where doublet labelling is stored
    OUT="${EXP_ID}/${SET}_filt_2ndRound"
    
    [[ -d "${OUT}" ]] || mkdir -p "${OUT}"

    # Run filtering by doublets and, optionally, by low-quality cluster
    # N.B.: the clustering have been computed on the original object, prior to doublet calculation, on the same dimensional reduction
    ID_FILT="filt_2ndRound_${EXP_ID}"
    COMMAND="${SCRIPT_DIR}/${FILT} \"${DIR_DOUBLETS}\" ${DR} \"${OUT}\" \"${DIR_CLUSTERS}\" ${SLOT} ${CL}" 
    sbatch -t 5:00:00 --mem=${MEM}G \
           -J ${ID_FILT} \
           -o ${OE}/${ID_FILT}.STDOUT \
           -e ${OE}/${ID_FILT}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"
}

# filter cells marked as doublets and, optionally, belonging to specified clusters

depend=""
for (( i=0; i< ${#EXP[@]}; i++ ))
do  
    exp=${EXP[$i]}
    if [[ ${CL_MODE[$i]} == "" ]]
    then
        # unsupervised: only remove doublets predicted by doubletFinder
        J=$( batch_job ${exp} ${exp}_ccRegress ${DR_DOUBL} )
        val=${J#Submitted batch job }
        depend="${depend}:${val}"
    else 
        # semi-supervised: remove a cluster of putative doublets, then remove the number of top-scoring doublets needed to reach the 10x estimate
        J=$( batch_job ${exp} ${exp}_ccRegress ${DR_DOUBL} ${CL_MODE[$i]} ${CL_FILT[$i]} )
        val=${J#Submitted batch job }
        depend="${depend}:${val}"
    fi

done
depend=${depend#:}

# generate the merged object and plot standard umap

ID="merged_N_filt_2ndRound_TIL"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/${MERGE} \"${DIR}\" ${IN_DIR_N[@]}"
sbatch -t 5:00:00 --mem=${MEM2}G \
       -J ${ID} \
       --dependency=afterok:${depend} \
       -o ${OE}/${ID}.STDOUT \
       -e ${OE}/${ID}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"

ID="merged_C_filt_2ndRound_no010"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/${MERGE} \"${DIR}\" ${IN_DIR_C_no010[@]}"
sbatch -t 5:00:00 --mem=${MEM2}G \
       -J ${ID} \
       --dependency=afterok:${depend} \
       -o ${OE}/${ID}.STDOUT \
       -e ${OE}/${ID}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"

ID="merged_all_filt_2ndRound"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/${MERGE} \"${DIR}\" ${IN_DIR[@]}"
sbatch -t 5:00:00 --mem=${MEM2}G \
       -J ${ID} \
       --dependency=afterok:${depend} \
       -o ${OE}/${ID}.STDOUT \
       -e ${OE}/${ID}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"


exit

