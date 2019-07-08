#$ -cwd
#job name
#$ -N STEP1_Haplotype_Concordance
#$ -o STEP1_Haplotype_Concordance.log
#$ -l mem=2G
#$ -l h_rt=10:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=2G
#$ -t 1-65

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

extension="bedfilter_no_indels_altref_altbam"
mapQfilter=30
minDP_S1=2
traw_file=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged.traw

Rscript Haplotype_Concordance.r ${extra_data_dir} ${traw_file} RData $SGE_TASK_ID


