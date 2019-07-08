#$ -cwd
#job name
#$ -N STEP1_prune_and_admixture
#$ -o STEP1_prune_and_admixture.log
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

#for plink, we speficy the allele to count so that we count alternate alleles:
awk -v OFS='\t' '{print $1":"$2,$5}' ${known_variants} | \
grep -v "^#" > ${TMPDIR}/reference_alleles.txt
#we will remove the outgroup and S2 samples from the genotypes used for ADMIXTURE analysis
remove_samples_file=${TMPDIR}/remove_samples_file.txt
echo "UC10164S2 UC10164S2" > $remove_samples_file
echo "Outgroup Outgroup" >> $remove_samples_file

for minDP_S1 in 2 4; do

plink_input=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged
plink_pruned_dir=${project_dir}/7Admixture/1LD_pruned_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/
plink_pruned=${plink_pruned_dir}/merged
admixture_out_dir=${project_dir}/7Admixture/2ADMIXTURE_out_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged
mkdir -p $plink_pruned_dir
mkdir -p $admixture_out_dir

plink \
--bfile $plink_input \
--remove $remove_samples_file \
--indep-pairwise 200 10 0.9 \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out $plink_pruned
plink \
--bfile $plink_input \
--remove $remove_samples_file \
--extract $plink_pruned.prune.in \
--make-bed \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${plink_pruned}.pruned

for rep in {1..50}; do 
  admixture_out_rep_dir=${admixture_out_dir}_rep${rep}
  mkdir -p $admixture_out_rep_dir
  cd $admixture_out_rep_dir
  
  for Kparam in {2..7}; do
    echo "admixture for K=$Kparam at $(date)"
    admixture --seed=1234${rep} ${plink_pruned}.pruned.bed $Kparam > merged.pruned.${Kparam}.rep${rep}.log 2>&1
    admixture --seed=1234${rep} --cv ${plink_pruned}.pruned.bed $Kparam | tee merged.pruned.${Kparam}.rep${rep}.CV.log
  done

done

done
