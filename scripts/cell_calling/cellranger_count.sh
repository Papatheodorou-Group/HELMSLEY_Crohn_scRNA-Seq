
# CR="$HOME/francesca-irenegroup/Software/cellranger-7.1.0"
# TSCRT="/nfs/research/marioni/common/references/hg38/refdata-gex-GRCh38-2020-A"

WDIR=$(pwd)
DATA_DIR=$( cd $(dirname $0) ; cd "../../data" )
TSCRT="${DATA_DIR}/references/hg38/refdata-gex-GRCh38-2020-A"
cd ${WDIR}

if [[ $# -lt "6" ]]
then
    echo "Usage: $0 <ncpu> <mem> <cell> <id> <samples> <fq>"
    exit
fi

NCPU=$1
MEM=$2
CELLS=$3
SAMPLE_NAME=$4
SAMPLE_ID=$5
FASTQS=$6

export PATH=$PATH:${CR}

cellranger count --nosecondary --localcores ${NCPU} --localmem ${MEM} --expect-cells ${CELLS} \
                 --id=${SAMPLE_NAME} --sample=${SAMPLE_ID} --transcriptome=${TSCRT} --fastqs=${FASTQS} --no-bam


exit

