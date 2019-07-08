#job name
#$ -N STEP1_Split_Raw_Fastq.reseq
#$ -o STEP1_Split_Raw_Fastq.reseq.log
# Set the working directory 
#$ -cwd
#memory and runtime options 
#$ -l mem=10G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=30G

SGE_TASK_ID=1

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#get all the fastqs, trim off extension that includes read pair, sort uniq to remove duplicate
for file in ${raw_fastq_dir}/*.fastq.gz ; do echo "${file%_R*.fastq.gz}"; done | sort -u > ${TMPDIR}/sample_list.txt

#select from sample_list the fastq file to work with
input_fastq_prefix=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/sample_list.txt)
#get the basename
experimentname=$(basename $input_fastq_prefix)
#trim the trailing part of the filename, including the lane and experiment number
samplename_in_experimentname=$(echo "${experimentname%_S[0-9]*}")

#the samplenames in the prefix are UC10164_S1 and UC10164S1. Convert both to get the name of the sample under consideration
if [[ $samplename_in_experimentname == *"S1"* ]]; then
  samplename=UC10164S1
else
  samplename=UC10164S2
fi

numthreads=1
num_lines_per_split=40000000
output_dir=${project_dir}/3AdapterRemoval/1Split_fastqs/${samplename}/${experimentname}
mkdir -p $output_dir

#split the fastq files using fastp 
fastp \
--in1 ${input_fastq_prefix}_R1_001.fastq.gz \
--in2 ${input_fastq_prefix}_R2_001.fastq.gz \
--out1 ${output_dir}/${experimentname}_R1.fastq.gz \
--out2 ${output_dir}/${experimentname}_R2.fastq.gz \
--thread $numthreads \
--split_by_lines $num_lines_per_split \
--disable_trim_poly_g \
--disable_adapter_trimming \
--disable_quality_filtering \
--disable_length_filtering \
--split_prefix_digits 4


