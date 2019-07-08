#$ -cwd
#job name
#$ -N STEP1_plot_Phylogeny_Dstats
#$ -o STEP1_plot_Phylogeny_Dstats.log
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

SNPhylo_dir=${project_dir}/6SNPhylo/2SNPhylo_result_${extension}/mapQ${mapQfilter}

for minDP_S1 in 2 4; do

for sub_dir in incl_UC10164S1 excl_UC10164S1; do

input_dir=${SNPhylo_dir}/S1minDP${minDP_S1}/${sub_dir}
Dstats_dir=${project_dir}/8AdmixTools/2AdmixTools_out_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/
plot_export_dir=plots/${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/${sub_dir}
mkdir -p $plot_export_dir

Rscript Phylogeny_Dstats.r $metadata_filename ${input_dir}/snphylo.output.bs.tree ${Dstats_dir} ${plot_export_dir}

done
done

sub_dir=excludeS1
input_dir=${SNPhylo_dir}/${sub_dir}
plot_export_dir=plots/${extension}/mapQ${mapQfilter}/${sub_dir}
mkdir -p $plot_export_dir

Rscript Phylogeny.r $metadata_filename ${input_dir}/snphylo.output.bs.tree ${plot_export_dir}

