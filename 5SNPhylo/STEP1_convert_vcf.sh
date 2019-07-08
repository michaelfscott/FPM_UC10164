#$ -cwd
#job name
#$ -N STEP1_convert_vcf
#$ -o STEP1_convert_vcf.log
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=20G
#$ -t 1-15

#parameters
max_depth_filter=30

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#get chromosome names and lengths for conversion
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
chr_names_lengths=${TMPDIR}/chromosome_names_lengths.txt
awk -v OFS='\t' '{print $1, $2, NR}' ${reference}.fai > ${chr_names_lengths}

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
merged_vcf=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/merged.vcf.gz

for minDP_S1 in 2 4; do

input_vcf=${TMPDIR}/merged_S1minDP${minDP_S1}.vcf
input_vcf2=${TMPDIR}/merged_S1minDP${minDP_S1}.excl_UC10164S1.vcf
outdir=${project_dir}/6SNPhylo/1converted_vcf_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}
mkdir -p ${outdir}
output_vcf=${outdir}/SNPhylo_incl_UC10164S1.vcf
output_vcf2=${outdir}/SNPhylo_excl_UC10164S1.vcf

#filter S1 depth and exclude UC10164S2
bcftools view $merged_vcf \
--include "FORMAT/DP[0] >= ${minDP_S1} && FORMAT/DP[0] < ${max_depth_filter}" \
-s ^UC10164S2 \
--output-type v \
--output-file $input_vcf
#filter S1 depth and exclude both UC10164
bcftools view $merged_vcf \
--include "FORMAT/DP[0] >= ${minDP_S1} && FORMAT/DP[0] < ${max_depth_filter}" \
--output-type v \
--output-file ${TMPDIR}/tmp.vcf
bcftools view ${TMPDIR}/tmp.vcf \
-s ^UC10164S1,UC10164S2 \
--output-type v \
--output-file $input_vcf2

#make a copy
cp $input_vcf ${output_vcf}.tmp.vcf

##### edit the names
while read line ; do
  chr=$(echo $line | awk '{print $1}')
  length=$(echo $line | awk '{print $2}')
  chrNum=$(echo $line | awk '{print $3}')
  echo "editing $chr of length $length at $(date)"
  awk -F $'\t' '$1 == "'$chr'" {$1 = "'$chrNum'"} {OFS = FS} {print}' ${output_vcf}.tmp.vcf > ${output_vcf}
  mv $output_vcf ${output_vcf}.tmp.vcf
done < ${chr_names_lengths}

mv ${output_vcf}.tmp.vcf ${output_vcf}

#bgzip and tabix index
bgzip --stdout $output_vcf > ${output_vcf}.gz
tabix ${output_vcf}.gz

#make a copy
cp $input_vcf2 ${output_vcf2}.tmp.vcf

##### edit the names
while read line ; do
  chr=$(echo $line | awk '{print $1}')
  length=$(echo $line | awk '{print $2}')
  chrNum=$(echo $line | awk '{print $3}')
  echo "editing $chr of length $length at $(date)"
  awk -F $'\t' '$1 == "'$chr'" {$1 = "'$chrNum'"} {OFS = FS} {print}' ${output_vcf2}.tmp.vcf > ${output_vcf2}
  mv $output_vcf2 ${output_vcf2}.tmp.vcf
done < ${chr_names_lengths}

mv ${output_vcf2}.tmp.vcf ${output_vcf2}

#bgzip and tabix index
bgzip --stdout $output_vcf2 > ${output_vcf2}.gz
tabix ${output_vcf2}.gz

done

#######
#now exclude UC10164S1 entirely
#######

input_vcf=${TMPDIR}/merged.vcf
outdir=${project_dir}/6SNPhylo/1converted_vcf_${extension}/mapQ${mapQfilter}/excludeS1
mkdir -p ${outdir}
output_vcf=${outdir}/SNPhylo_excl_UC10164S1.vcf

#filter S1 and S2
bcftools view $merged_vcf \
-s ^UC10164S1,UC10164S2 \
--output-type v \
--output-file $input_vcf

#make a copy
cp $input_vcf ${output_vcf}.tmp.vcf

##### edit the names
while read line ; do
  chr=$(echo $line | awk '{print $1}')
  length=$(echo $line | awk '{print $2}')
  chrNum=$(echo $line | awk '{print $3}')
  echo "editing $chr of length $length at $(date)"
  awk -F $'\t' '$1 == "'$chr'" {$1 = "'$chrNum'"} {OFS = FS} {print}' ${output_vcf}.tmp.vcf > ${output_vcf}
  mv $output_vcf ${output_vcf}.tmp.vcf
done < ${chr_names_lengths}

mv ${output_vcf}.tmp.vcf ${output_vcf}

#bgzip and tabix index
bgzip --stdout $output_vcf > ${output_vcf}.gz
tabix ${output_vcf}.gz

