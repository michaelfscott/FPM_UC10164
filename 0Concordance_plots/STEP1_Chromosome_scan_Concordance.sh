#$ -cwd
#job name
#$ -N STEP1_Chromosome_scan_Concordance
#$ -o STEP1_Chromosome_scan_Concordance.log
#$ -l mem=10G
#$ -l h_rt=10:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=2G

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

extension="bedfilter_no_indels_altref_altbam"
mapQfilter=30
minDP_S1=2
traw_file=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged.traw
traw_regions_file=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/merged.regions.traw

Rscript Chromosome_scan_Concordance.r ${extra_data_dir} ${traw_file} ${traw_regions_file} plots/


