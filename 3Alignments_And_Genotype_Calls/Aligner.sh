#parameters
numthreads=4
bwa_seedlength=16500
bwa_editdist=0.01
bwa_gapopens=2

#we look in the datadir for the single end reads after adapter removal
ls ${datadir}/*/*/*.collapsed.not_truncated.fastq.gz > ${TMPDIR}/SE_${read_type}_fastq_list.txt

#select from sample_list the fastq file to work with
fastq=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/SE_${read_type}_fastq_list.txt)
prefix=$( basename $fastq .collapsed.not_truncated.fastq.gz )

#The fastq files should be in subdirectories named by the samplename and the experiment name, this extracts those parts of the path. 
samplename=$(echo $fastq | awk -F'/' '{print $(NF-2)}')
experimentname=$(echo $fastq | awk -F'/' '{print $(NF-1)}')

#create output directory
outdir=${project_dir}/4Alignments/1Align_${species}_${read_type}/SplitReads/${samplename}/
mkdir -p ${outdir}

#extract the read group information from the fastq file
string=$(gzip -cd $fastq | head -n 2 | grep "^@") #get the read name from the head of fastq1
instrument_tmp=$(echo $string | cut -d":" -f 1) #first field is a unique instrument name
instrument=${instrument_tmp#"@"} #trim off leading @ from the fastq read string
flowcell=$(echo $string | cut -d":" -f 3) #flowcell ID
lane=$(echo $string | cut -d":" -f 4) #the lane is the fourth part.
#the barcode is different for UC10164S1 and UC10164S2, manually entered here
if [[ $samplename == *"UC10164S1"* ]]; then
  barcode=CAATTAC
else
  barcode=AGATAGG
fi
#there is only one library sequenced for each sample
library=LIB${samplename}
#put this information together to get the RG for bwa
runRG="@RG\tID:${flowcell}.${lane}\tLB:${library}\tPL:ILLUMINA\tSM:${samplename}\tPU:${flowcell}.${lane}.${barcode}"

#make a temporary directory for samtools sort
mkdir -p ${TMPDIR}/sort_${prefix}

echo "$SGE_TASK_ID aligning SE reads $prefix at $(date)"
bwa aln -t $numthreads -l $bwa_seedlength -n $bwa_editdist -o $bwa_gapopens $reference ${fastq} > ${TMPDIR}/${prefix}_sa.sai 
echo "$SGE_TASK_ID samse and samtools sort for $prefix at $(date)"
bwa samse -r ${runRG} $reference ${TMPDIR}/${prefix}_sa.sai ${fastq} | \
samtools sort -@ $numthreads -T ${TMPDIR}/sort_${prefix}/${prefix} -o ${TMPDIR}/${prefix}.sorted.bam -
echo "$SGE_TASK_ID indexing ${prefix} at $(date)"
samtools index ${TMPDIR}/${prefix}.sorted.bam

echo "$SGE_TASK_ID copying at $(date)"
cp ${TMPDIR}/${prefix}.sorted.bam* --target-directory=${outdir}/

