source ~/.bashrc
conda activate scvi-env

if [[ $# -lt "3" ]]
then
    echo "usage: bash plot_celltypes.sh <in_obj> <out_dir> <reduct>"
    exit
fi

SCRIPT_DIR=$( cd $( dirname $0 ) ; pwd )
IN_OBJ=$1
OUT_DIR=$2
REDUCT=$3

python ${SCRIPT_DIR}/plot_celltypes.py "${IN_OBJ}" "${OUT_DIR}" ${REDUCT}

exit

