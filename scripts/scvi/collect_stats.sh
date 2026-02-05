
META=( "category" "Integrated_05" ) # broad or fine-grained cell-type annotation

if [[ $# -lt "2" ]]
then
    echo "usage: bash collect_stats.sh <r_prefix> <out>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
R_OBJ_PREFIX=$1
OUT=$2

Rscript ${SCRIPT_DIR}/collect_stats.R "${R_OBJ_PREFIX}" "${OUT}" ${META[@]}


exit

