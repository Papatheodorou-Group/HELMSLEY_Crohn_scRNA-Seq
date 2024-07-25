
# module purge
# module load r-4.1.0-gcc-9.3.0-wvnko7v gmp-6.1.2-gcc-9.3.0-hicntdj
# R_LIBS_USER="/hps/software/users/marioni/francesca/R_libs"
# export R_LIBS_USER

SCRIPT_DIR=$( cd $(dirname $0) ; pwd )

if [[ $# -lt 4 ]]
then
    echo "Usage: bash $0 <dir> <mode> <cl> <outdir>"
    exit
fi

dir=$1
mode=$2
cl=$3
outdir=$4

[[ -d ${outdir} ]] || mkdir -p ${outdir}
cells="${outdir}/cells.txt"

# extract the cells
Rscript ${SCRIPT_DIR}/cell_subset.R ${dir}/hvg_pca_clust/object.Rds ${mode} "${cl}" ${cells}

# subset the object
Rscript ${SCRIPT_DIR}/obj_subset.R ${dir}/object.Rds ${cells} ${outdir}/object.Rds

exit

