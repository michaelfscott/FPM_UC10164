#$ -cwd
#job name
#$ -N STEP1_Align_Emmer
#$ -o STEP1_Align_Emmer.log
#memory and runtime options 
#$ -l mem=8G
#Number of threads per job, memory given above is per-thread
#$ -pe smp 4
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=60G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-87

#export variables
set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
species=Emmer
read_type=Collapsed
datadir=${project_dir}/3AdapterRemoval/2Adapters_Removed

bash Aligner.sh
