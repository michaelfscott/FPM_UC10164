#job name
#$ -N STEP1_Split_Bread_Reference
#$ -o STEP1_Split_Bread_Reference.log
# Set the working directory 
#$ -cwd
#memory and runtime options 
#$ -l mem=40G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=30G

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

java_memory_param=10

input_reference_genome=${bread_wheat_reference_original}
output_dir=${project_dir}/2Reference_Genomes/IWGSC_refseq_v1_split
output_file=${output_dir}/161010_Chinese_Spring_v1.0_pseudomolecules_split.fasta
split_fasta_dir=${TMPDIR}/separate_chromosomes

mkdir -p $output_dir
mkdir -p $split_fasta_dir

#index the input reference genome
echo "samtools faidx emmer reference at $(date)"
samtools faidx $input_reference_genome

#The fai file has the names and lengths of all the chromosomes
awk '{print $1 "\t" $2}' ${input_reference_genome}.fai > ${TMPDIR}/chr_names_lengths.txt

#for each chromosome, split into two halves. ${chr}_1 has the positions from 1-400mb, ${chr}_2 has the remaining positions.
while read line ; do
  chr=$(echo $line | awk '{print $1}')
  length=$(echo $line | awk '{print $2}')
  echo $chr
  echo $length
  samtools faidx $input_reference_genome ${chr}:1-400000000 > ${split_fasta_dir}/${chr}_1.fasta 
  samtools faidx $input_reference_genome ${chr}:400000001-${length} > ${split_fasta_dir}/${chr}_2.fasta 
done < ${TMPDIR}/chr_names_lengths.txt

#join all split chromosome halves into one fasta file to get new split reference genome
cat ${split_fasta_dir}/*.fasta > ${output_file}

#index reference genome
echo "bwa indexing reference at $(date)"
bwa index ${output_file}
echo "samtools faidx reference at $(date)"
samtools faidx ${output_file}
echo "GATK indexing reference at $(date)"
gatk CreateSequenceDictionary --java-options "-Xmx${java_memory_param}G" \
-R ${output_file}

