#$ -cwd
#job name
#$ -N STEP2_Merge_Emmer_Trimmed
#$ -o STEP2_Merge_Emmer_Trimmed.log
#memory and runtime options 
#$ -l mem=30G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=150G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-2

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#parameters
numthreads=1
java_memory_param=18

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
read_type=Collapsed_Trimmed

bash Merge_Bams_copy_fullbam.sh

