#job name
#$ -N STEP2_AdapterRemoval_And_Trimming.blank
#$ -o STEP2_AdapterRemoval_And_Trimming.blank.log
# Set the working directory 
#$ -cwd
#memory and runtime options 
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=10G
#request number of threads (each has different SGE_TASK_ID variable)

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

samplename=blank
experimentname=Sel_II_S2

#the split paired end reads should be in the following directory (from STEP1)
ls ${project_dir}/1raw_Blank_fastq/Sel_II/*_R1_001.fastq.gz > ${TMPDIR}/R1_fastq_list.txt
ls ${project_dir}/1raw_Blank_fastq/Sel_II/*_R2_001.fastq.gz > ${TMPDIR}/R2_fastq_list.txt

for SGE_TASK_ID in {1..4}; do

#select from sample_list the fastq file to work with
fastq_R1=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/R1_fastq_list.txt)
fastq_R2=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/R2_fastq_list.txt)

output_dir=${project_dir}/3AdapterRemoval/2Adapters_Removed/${samplename}/${experimentname}
output_dir_trimmed=${project_dir}/3AdapterRemoval/3Adapters_Removed_And_Trimmed/${samplename}/${experimentname}
mkdir -p $output_dir
mkdir -p ${output_dir_trimmed}

prefix=$(basename $fastq_R1 _R1_001.fastq.gz)
output_base=${output_dir}/${prefix}

numthreads=1

echo "adapter removal for $prefix begins at $(date)"

AdapterRemoval --file1 ${fastq_R1} --file2 ${fastq_R2} \
--threads ${numthreads} \
--gzip \
--minlength 20 \
--trimns \
--trimqualities \
--collapse \
--adapter1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT \
--basename ${output_base} \
--output1 ${output_base}.mate1.fastq.gz \
--output2 ${output_base}.mate2.fastq.gz \
--outputcollapsed ${output_base}.collapsed.not_truncated.fastq.gz \
--discarded ${output_base}.discarded.fastq.gz \
--singleton ${output_base}.singleton.fastq.gz \
--settings ${output_base}.adapterremoval.settings 

output_base_trimmed=${output_dir_trimmed}/${prefix}

echo trimming $prefix at $(date)
zcat ${output_base}.collapsed.not_truncated.fastq.gz | \
fastx_trimmer -Q33 -f 3 | \
fastx_trimmer -Q33 -t 2 | \
gzip > ${output_base_trimmed}.trimmed.collapsed.not_truncated.fastq.gz

done
