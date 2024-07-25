
# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

if [[ $# -lt "5" ]]
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

OBJ="${DIR}/object.Rds"

META_COL=$( echo ${META_LIST[@]} | tr " " "," )

# add metadata and update the object
Rscript ${GIT_DIR}/Seurat_add_meta_data.R --object ${OBJ} --table ${METAFILE} --columns ${META_COL}


exit
