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

SUFFIX="ccRegress_filt"
IN_DIR_N=( $( for exp in ${EXP_N[@]} ; do echo ${exp}/${exp}_${SUFFIX} ; done ) )
IN_DIR_C=( $( for exp in ${EXP_C[@]} ; do echo ${exp}/${exp}_${SUFFIX} ; done ) )
IN_DIR_C_no010=( $( for exp in ${EXP_C_no010[@]} ; do echo ${exp}/${exp}_${SUFFIX} ; done ) )
IN_DIR=( ${IN_DIR_N[@]} ${IN_DIR_C_no010[@]} )

# select the clustering solution with at least 20 clusters and best silhouette score (vst5000, k=30 only)
# check mitochondrial gene expression in clusters from DEA: if it is not clear-cut, increase clustering resolution
CL_MODE=( "clusters_pca_vst_top5000_k30_res0.6" \
          "clusters_pca_vst_top5000_k30_res0.8" \
          "clusters_pca_vst_top5000_k30_res0.8" \
          "clusters_pca_vst_top5000_k30_res1.2" \
          "clusters_pca_vst_top5000_k30_res1" \
          "clusters_pca_vst_top5000_k30_res0.8" \
          "clusters_pca_vst_top5000_k30_res0.7" \
          "clusters_pca_vst_top5000_k30_res1.2" ) 

# clusters to remove
CL_FILT=( "11" "11" "4" "2,8" "14" "10" "4" "" ) 

# as per 10x v3 guidelines
DOUBLET_PERC=( "7.6" "7.6" "7.6" "7.6" "7.6" "10" "7.6" "7.6" ) # assume ~10000 expected cells for each sample
                                                                # increased the doublet rate for 008 due to 10x machine overloading
       
MEM="30"
MEM2="50"

FILT="cell_filtering.sh"
DOUBL="find_doublets.sh"
MERGE="merge.sh"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}

# filter cells belonging to specified clusters

depend=""
for (( i=0; i<${#EXP[@]}; i++ ))
do  
    EXP_ID=${EXP[$i]}
    SET="${EXP[$i]}_ccRegress"
    SLOT=${CL_MODE[$i]}
    CL=${CL_FILT[$i]}
    DP=${DOUBLET_PERC[$i]}

    DIR="${EXP_ID}/${SET}"
    OUT="${EXP_ID}/${SET}_filt"
    
    [[ -d "${OUT}" ]] || mkdir -p "${OUT}"
    
    # Run filtering by low-quality cluster
    ID_FILT="filt_${EXP_ID}"
    COMMAND="${SCRIPT_DIR}/${FILT} \"${DIR}\" ${SLOT} \"${CL}\" \"${OUT}\"" 
    J=$(sbatch -t 5:00:00 --mem=${MEM}G \
               -J ${ID_FILT} \
               -o ${OE}/${ID_FILT}.STDOUT \
               -e ${OE}/${ID_FILT}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")
    J1=${J#Submitted batch job }
    
    # Re-process and find doublets
    # save the doublet annotation to a tsv file
    # collect the count/perc of doublets in clusters 
    ID_DOUBL="doubl_${EXP_ID}"
    COMMAND="${SCRIPT_DIR}/${DOUBL} \"${GIT_DIR}\" \"${DIR}\" \"${OUT}\" \"${WDIR}\" ${DP} ${SLOT}"
    J=$(sbatch -t 5:00:00 --mem=${MEM}G \
               -J ${ID_DOUBL} \
               --dependency=afterok:${J1} \
               -o ${OE}/${ID_DOUBL}.STDOUT \
               -e ${OE}/${ID_DOUBL}.STDERR \
               --wrap="${PREAMBLE} ; bash ${COMMAND}")              
    J2=${J#Submitted batch job }  

    depend="${depend}:${J2}"
done
depend=${depend#:}

# generate the merged objects and plot standard umap
# N.B.: merged objects do not contain cluster nor doublet annotation of the individual samples

ID="merged_N_filt_TIL"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/${MERGE} \"${DIR}\" ${IN_DIR_N[@]}"
sbatch -t 5:00:00 --mem=${MEM2}G \
       -J ${ID} \
       --dependency=afterok:${depend} \
       -o ${OE}/${ID}.STDOUT \
       -e ${OE}/${ID}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"

ID="merged_C_filt_no010"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/${MERGE} \"${DIR}\" ${IN_DIR_C_no010[@]}"
sbatch -t 5:00:00 --mem=${MEM2}G \
       -J ${ID} \
       --dependency=afterok:${depend} \
       -o ${OE}/${ID}.STDOUT \
       -e ${OE}/${ID}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"
              
ID="merged_all_filt"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/${MERGE} \"${DIR}\" ${IN_DIR[@]}"
sbatch -t 5:00:00 --mem=${MEM2}G \
       -J ${ID} \
       --dependency=afterok:${depend} \
       -o ${OE}/${ID}.STDOUT \
       -e ${OE}/${ID}.STDERR \
       --wrap="${PREAMBLE} ; bash ${COMMAND}"

exit

