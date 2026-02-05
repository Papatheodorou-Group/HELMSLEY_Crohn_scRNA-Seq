# run cellranger 
# cell recovery based on https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf

LOCALDIR=$( cd $( dirname $0 ) ; pwd )
SCRIPTDIR=$( cd "../../scripts/cell_calling/" ; pwd )

NCPU="4"
MEM="64"

SAMPLE_INFO="sample_info.tsv"

IFS=$'\n'
for line in $(cat ${SAMPLE_INFO})
do
    SAMPLE_NAME=$(echo $line | cut -f 1)
    SAMPLE_ID=$(echo $line | cut -f 2)
    FASTQS=${LOCALDIR}/$(echo $line | cut -f 3)
    CELLS=$(echo $line | cut -f 4) 

    ID="cr_count_${SAMPLE_NAME}"
    COMMAND="${SCRIPTDIR}/cellranger_count.sh ${NCPU} ${MEM} ${CELLS} ${SAMPLE_NAME} ${SAMPLE_ID} ${FASTQS}"
    sbatch -t 48:00:00 --mem=${MEM}G -c ${NCPU} \
           -J ${ID} \
           -o ${ID}.STDOUT \
           -e ${ID}.STDERR \
           --wrap="bash ${COMMAND}"
done


exit



