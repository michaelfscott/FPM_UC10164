#$ -cwd
#job name
#$ -N STEP2_get_max_LogLik
#$ -o STEP2_get_max_LogLik.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=20G
#$ -t 1-15

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

#choose mapQ filter and samplename parameters across SGE_TASK_ID's
get_filters() {
  local TASK_ID=$1
  local mapQs="20 25 30"
  local filters=".bedfilter .bedfilter.noindels .bedfilter.posfilter.altref .bedfilter.noindels.posfilter.altBam .bedfilter.noindels.posfilter.altRef.altBam"
  local filter_extensions="bedfilter bedfilter_no_indels bedfilter_altref bedfilter_no_indels_altbam bedfilter_no_indels_altref_altbam"
  local mapQ_num=$((($TASK_ID+2)%3 + 1))
  mapQfilter=$(echo $mapQs | cut -d" " -f$mapQ_num)
  local filter_num=$((($TASK_ID+2)/3))
  filter=$(echo $filters | cut -d" " -f$filter_num)
  extension=$(echo $filter_extensions | cut -d" " -f$filter_num)
}

get_filters $SGE_TASK_ID

for minDP_S1 in 2 4; do

admixture_out_dir=${project_dir}/7Admixture/2ADMIXTURE_out_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged
mkdir -p $admixture_out_dir
rm ${admixture_out_dir}/merged.pruned*

for K in {2..7}; do
	rep_CV_filename=$(grep -A 2 "Summary" ${admixture_out_dir}*/merged.pruned.${K}.rep*.CV.log | grep Loglik | sort -nrk2,2 | head -n 1 | cut -f1 -d" " | cut -f1 -d"-")
	rep_prefix=${rep_CV_filename%*rep*}
	cp ${rep_prefix}Q $admixture_out_dir
	cp ${rep_CV_filename} ${admixture_out_dir}
done
grep "CV" ${admixture_out_dir}/*.CV.log | awk '{print $NF}' > ${admixture_out_dir}/CV_error.txt

done
