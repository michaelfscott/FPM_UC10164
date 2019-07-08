#$ -cwd
#job name
#$ -N STEP2_Haplotype_Concordance
#$ -o STEP2_Haplotype_Concordance.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=2G
#$ -t 1-10

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

extension="bedfilter_no_indels_altref_altbam"
mapQfilter=30
minDP_S1=2
traw_file=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged.traw

Rscript Haplotype_Concordance.2.r ${extra_data_dir} ${traw_file} $SGE_TASK_ID

