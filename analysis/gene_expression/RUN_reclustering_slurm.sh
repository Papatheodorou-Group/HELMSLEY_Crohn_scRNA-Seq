
#### PARAMETERS #####

PREAMBLE="module purge
module load r-4.2.2-gcc-11.2.0-oa3uudy gmp-6.2.1-gcc-11.2.0-mneucsf 
R_LIBS_USER=/hps/software/users/marioni/francesca/R_libs
export R_LIBS_USER"

SCRIPT_DIR=$( cd "../../scripts/seurat_analysis" ; pwd )
SCANPY_DIR=$( cd "../annotation/HELMSLEY/stats" ; pwd )
GCA_DIR=$( cd "../annotation/GCA" ; pwd )
GIT_DIR="../../../r_scripts/gene_expression_analysis"
WDIR=$(pwd)

EXP_N=( "004_N_E_S" \
        "006_N_E_S" \
        "007_N_E_S" \
        "012_N_E_S" )
EXP_C=( "003_C_E_S" \
        "008_C_E_S" \
        "013_C_E_S" )
EXP=( ${EXP_N[@]} ${EXP_C[@]} )      

SUFFIX_ALL="ccRegress_filt"
IN_DIR_N=( $( for exp in ${EXP_N[@]} ; do echo ${exp}/${exp}_${SUFFIX_ALL} ; done ) )
IN_DIR_C=( $( for exp in ${EXP_C[@]} ; do echo ${exp}/${exp}_${SUFFIX_ALL} ; done ) )
IN_DIR=( ${IN_DIR_N[@]} ${IN_DIR_C[@]} )

SUFFIX_ALL_2nd="ccRegress_filt_2ndRound"
IN_DIR_N_2nd=( $( for exp in ${EXP_N[@]} ; do echo ${exp}/${exp}_${SUFFIX_ALL_2nd} ; done ) )
IN_DIR_C_2nd=( $( for exp in ${EXP_C[@]} ; do echo ${exp}/${exp}_${SUFFIX_ALL_2nd} ; done ) )
IN_DIR_2nd=( ${IN_DIR_N[@]} ${IN_DIR_C[@]} )

GCA_STATS="${GCA_DIR}/GCA_stats.tsv"

META_SCANPY=( "category" "Integrated_05" "category_inferred" )
META_DOUBLET=( "doublet_label" "doublet_score" )
SCANPY_INFO="${SCANPY_DIR}/merged_all_filt.tsv"

DOUBLET_INFO_N=( $( for dir in ${IN_DIR_N[@]} ; do echo "${dir}/doublet_res.tsv" ; done ) )
DOUBLET_INFO_C=( $( for dir in ${IN_DIR_C[@]} ; do echo "${dir}/doublet_res.tsv" ; done ) )
                         
MEM="60"

OE="oe"
[[ -d $OE ]] || mkdir ${OE}



#### FUNCTIONS ####

batch_job_1 () {

    SET=$1
    LABEL=$2
    SUFFIX=$3
    METAFILE1=$4
    METAFILE2=$5
    DEPEND=$6

    OUT_OBJ="${LABEL}/${LABEL}_${SUFFIX}"
    [[ -d "${OUT_OBJ}" ]] || mkdir -p "${OUT_OBJ}"

    # Add scANVI metadata to the Seurat object
    ID1="add_meta_${LABEL}"
    COMMAND="${SCRIPT_DIR}/add_metadata_scanvi.sh \"${GIT_DIR}\" \"${METAFILE1}\" \"${SET}\" \"${OUT_OBJ}\" ${META_SCANPY[@]}"
    if [[ ${DEPEND} == "" ]]
    then
        J=$(sbatch -t 2:00:00 --mem=${MEM}G \
                   -J ${ID1} \
                   -o ${OE}/${ID1}.STDOUT \
                   -e ${OE}/${ID1}.STDERR \
                   --wrap="${PREAMBLE} ; bash ${COMMAND}")
    else
        J=$(sbatch -t 2:00:00 --mem=${MEM}G \
                   -J ${ID1} \
                   --dependency=afterok:${DEPEND} \
                   -o ${OE}/${ID1}.STDOUT \
                   -e ${OE}/${ID1}.STDERR \
                   --wrap="${PREAMBLE} ; bash ${COMMAND}")
    fi
    J1=${J#Submitted batch job }
    
    # Add doublet information
    ID2="add_doublet_${LABEL}"
    COMMAND="${SCRIPT_DIR}/add_metadata_df.sh \"${GIT_DIR}\" \"${METAFILE2}\" \"${SET}\" \"${OUT_OBJ}\" ${META_DOUBLET[@]}"
    sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID2} \
           --dependency=afterok:${J1} \
           -o ${OE}/${ID2}.STDOUT \
           -e ${OE}/${ID2}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"
}

batch_job_2 () {

    SET=$1
    LABEL=$2
    SUFFIX=$3
    DEPEND=$4
    
    OUT_OBJ="${LABEL}/${LABEL}_${SUFFIX}"
    [[ -d "${OUT_OBJ}" ]] || mkdir -p "${OUT_OBJ}"

    # Run the clustering
    ID3="reclust_${LABEL}"
    COMMAND="${SCRIPT_DIR}/5_reclustering.sh \"${GIT_DIR}\" \"${OUT_OBJ}\" \"${WDIR}\"" 
    if [[ ${DEPEND} == "" ]]
    then
        J=$(sbatch -t 10:00:00 --mem=${MEM}G \
                   -J ${ID3} \
                   -o ${OE}/${ID3}.STDOUT \
                   -e ${OE}/${ID3}.STDERR \
                   --wrap="${PREAMBLE} ; bash ${COMMAND}")
    else
        J=$(sbatch -t 10:00:00 --mem=${MEM}G \
                   -J ${ID3} \
                   --dependency=afterok:${DEPEND} \
                   -o ${OE}/${ID3}.STDOUT \
                   -e ${OE}/${ID3}.STDERR \
                   --wrap="${PREAMBLE} ; bash ${COMMAND}")
    fi
    J3=${J#Submitted batch job }
   
    # collect statistics on generated clustering solutions
    ID4="stats_${LABEL}"
    COMMAND="${SCRIPT_DIR}/6_collect_recl_stat.sh \"${WDIR}\" \"${OUT_OBJ}\" \"${GCA_STATS}\"" 
    sbatch -t 2:00:00 --mem=${MEM}G \
           -J ${ID4} \
           --dependency=afterok:${J3} \
           -o ${OE}/${ID4}.STDOUT \
           -e ${OE}/${ID4}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}"          
}



### EXECUTE ####


## without doublet filtering

dependN=""
for (( i=0; i<${#EXP_N[@]}; i++ ))
do
    J=$( batch_job_1 ${EXP_N[$i]} ${EXP_N[$i]} ${SUFFIX_ALL} ${SCANPY_INFO} ${DOUBLET_INFO_N[$i]} )
    val1=${J#Submitted batch job }
    J=$( batch_job_2 ${EXP_N[$i]} ${EXP_N[$i]} ${SUFFIX_ALL} ${val1} )
    val2=${J#Submitted batch job }
    dependN="${dependN}:${val2}"
done
dependN=${dependN#:}

dependC=""
for (( i=0; i<${#EXP_C[@]}; i++ ))
do
    J=$( batch_job_1 ${EXP_C[$i]} ${EXP_C[$i]} ${SUFFIX_ALL} ${SCANPY_INFO} ${DOUBLET_INFO_C[$i]} )
    val1=${J#Submitted batch job }
    J=$( batch_job_2 ${EXP_C[$i]} ${EXP_C[$i]} ${SUFFIX_ALL} ${val1} )
    val2=${J#Submitted batch job }
    dependC="${dependC}:${val2}"
done
dependC=${dependC#:}

# merge normal TIL samples containing doublet information
echo "merge normal TIL samples containing doublet information"
ID="merged_N_filt_TIL"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/merge.sh \"${DIR}\" ${IN_DIR_N[@]}"
J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           --dependency=afterok:${dependN} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J4=${J#Submitted batch job }

# cluster the merged samples
N_SAMPLES=$( echo ${EXP_N[@]} | tr ' ' ',' )
batch_job_2 ${N_SAMPLES} "merged" "N_filt_TIL" ${J4}

# merge Crohn samples containing doublet information
echo "merge Crohn samples containing doublet information"
ID="merged_C_filt_no010"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/merge.sh \"${DIR}\" ${IN_DIR_C[@]}"
J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           --dependency=afterok:${dependC} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J5=${J#Submitted batch job }

# cluster the merged samples
C_SAMPLES=$( echo ${EXP_C[@]} | tr ' ' ',' )
batch_job_2 ${C_SAMPLES} "merged" "C_filt_no010" ${J5}

# merge all samples containing doublet information
echo "merge all samples containing doublet information"
ID="merged_all_filt"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/merge.sh \"${DIR}\" ${IN_DIR[@]}"
J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           --dependency=afterok:${dependN}:${dependC} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J6=${J#Submitted batch job }

# cluster the merged samples
ALL_SAMPLES=$( echo ${EXP[@]} | tr ' ' ',' )
batch_job_2 ${ALL_SAMPLES} "merged" "all_filt" ${J6}




## with doublet filtering

dependN=""
for (( i=0; i<${#EXP_N[@]}; i++ ))
do
    J=$( batch_job_1 ${EXP_N[$i]} ${EXP_N[$i]} ${SUFFIX_ALL_2nd} ${SCANPY_INFO} ${DOUBLET_INFO_N[$i]} )
    val1=${J#Submitted batch job }
    J=$( batch_job_2 ${EXP_N[$i]} ${EXP_N[$i]} ${SUFFIX_ALL_2nd} ${val1} )
    val2=${J#Submitted batch job }
    dependN="${dependN}:${val2}"
done
dependN=${dependN#:}

dependC=""
for (( i=0; i<${#EXP_C[@]}; i++ ))
do
    J=$( batch_job_1 ${EXP_C[$i]} ${EXP_C[$i]} ${SUFFIX_ALL_2nd} ${SCANPY_INFO} ${DOUBLET_INFO_C[$i]} )
    val1=${J#Submitted batch job }
    J=$( batch_job_2 ${EXP_C[$i]} ${EXP_C[$i]} ${SUFFIX_ALL_2nd} ${val1} )
    val2=${J#Submitted batch job }
    dependC="${dependC}:${val2}"
done
dependC=${dependC#:}

# merge normal TIL samples containing doublet information
echo "merge normal TIL samples with doublets filtered out"
ID="merged_N_filt_2ndRound_TIL"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/merge.sh \"${DIR}\" ${IN_DIR_N_2nd[@]}"
J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           --dependency=afterok:${dependN} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J4=${J#Submitted batch job }

# cluster the merged samples
N_SAMPLES=$( echo ${EXP_N[@]} | tr ' ' ',' )
batch_job_2 ${N_SAMPLES} "merged" "N_filt_2ndRound_TIL" ${J4}

# merge Crohn samples containing doublet information
echo "merge Crohn samples with doublets filtered out"
ID="merged_C_filt_2ndRound_no010"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/merge.sh \"${DIR}\" ${IN_DIR_C_2nd[@]}"
J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           --dependency=afterok:${dependC} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J5=${J#Submitted batch job }

# cluster the merged samples
C_SAMPLES=$( echo ${EXP_C[@]} | tr ' ' ',' )
batch_job_2 ${C_SAMPLES} "merged" "C_filt_2ndRound_no010" ${J5}

# merge all samples containing doublet information
echo "merge all samples with doublets filtered out"
ID="merged_all_filt_2ndRound"
DIR="merged/${ID}"
COMMAND="${SCRIPT_DIR}/merge.sh \"${DIR}\" ${IN_DIR_2nd[@]}"
J=$(sbatch -t 12:00:00 --mem=${MEM}G \
           -J ${ID} \
           --dependency=afterok:${dependN}:${dependC} \
           -o ${OE}/${ID}.STDOUT \
           -e ${OE}/${ID}.STDERR \
           --wrap="${PREAMBLE} ; bash ${COMMAND}")
J6=${J#Submitted batch job }     

# cluster the merged samples
ALL_SAMPLES=$( echo ${EXP[@]} | tr ' ' ',' )
batch_job_2 ${ALL_SAMPLES} "merged" "all_filt_2ndRound" ${J6}


exit

