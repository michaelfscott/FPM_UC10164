#$ -cwd
#job name
#$ -N STEP1_plot_ADMIXTURE
#$ -o STEP1_plot_ADMIXTURE.log
#$ -l mem=2G
#$ -l h_rt=2:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=2G
#$ -t 1-3

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

extension="bedfilter_no_indels_altref_altbam"
mapQ_num=$SGE_TASK_ID
mapQs="20 25 30"
mapQfilter=$(echo $mapQs | cut -d" " -f$mapQ_num)

metadata_filename=${extra_data_dir}/origins_with_location.txt

for minDP_S1 in 2 4; do

plink_pruned_dir=${project_dir}/7Admixture/1LD_pruned_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/
plink_pruned=${plink_pruned_dir}/merged
admixture_out_dir=${project_dir}/7Admixture/2ADMIXTURE_out_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged
plot_export_dir=plots/${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}
mkdir -p $plot_export_dir

Rscript plot_ADMIXTURE.r $metadata_filename ${plink_pruned}.pruned.fam ${admixture_out_dir} ${plot_export_dir}

done
