#$ -cwd
#job name
#$ -N STEP2_Merge_Bread
#$ -o STEP2_Merge_Bread.log
#memory and runtime options 
#$ -l mem=30G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=150G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-2

#export variables
set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#parameters
numthreads=1
java_memory_param=18

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/IWGSC_refseq_v1_split/161010_Chinese_Spring_v1.0_pseudomolecules_split.fasta
species=Bread
read_type=Collapsed

bash Merge_Bams_copy_fullbam.sh

