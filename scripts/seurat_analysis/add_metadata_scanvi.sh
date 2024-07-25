
# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

if [[ $# -lt "6" ]]
then
    echo ""
    echo "Usage: bash $0 <git_dir> <metafile> <sample_name> <dir> <meta1> ... <metaN>"
    echo ""
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )

GIT_DIR=$1
METAFILE=$2
SAMPLE_NAME=$3
DIR=$4

n=$(($#-4))
META_LIST=( ${@:5:$n} )

META="${DIR}/metadata_scanvi.tsv"
OBJ="${DIR}/object.Rds"

META_COL=$( echo ${META_LIST[@]} | tr " " "," )

# reformat the table
Rscript ${SCRIPT_DIR}/reformat_meta_from_scanvi.R ${METAFILE} ${SAMPLE_NAME} ${META} ${META_COL}

# add metadata and update the object
Rscript ${GIT_DIR}/Seurat_add_meta_data.R --object ${OBJ} --table ${META} --columns ${META_COL}


exit
