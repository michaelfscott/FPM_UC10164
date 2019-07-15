#$ -cwd
#job name
#$ -N STEP7_Merge_VCFs_all
#$ -o STEP7_Merge_VCFs_all.log
#memory and runtime options 
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=20G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-15

#parameters
numthreads=1
max_depth_filter=30

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#choose mapQ filter and samplename parameters across SGE_TASK_ID's
get_filters() {
  local TASK_ID=$1
  #in order to cycle through the chromosomes and subgenomes
  local mapQs="20 25 30"
#  local samplenames="UC10164S1 UC10164S2"
  local filters=".bedfilter .bedfilter.noindels .bedfilter.posfilter.altref .bedfilter.noindels.posfilter.altBam .bedfilter.noindels.posfilter.altRef.altBam"
  local filter_extensions="bedfilter bedfilter_no_indels bedfilter_altref bedfilter_no_indels_altbam bedfilter_no_indels_altref_altbam"
  local mapQ_num=$((($TASK_ID+2)%3 + 1))
  mapQfilter=$(echo $mapQs | cut -d" " -f$mapQ_num)
#  local samplename_tmp=$((($TASK_ID +2)/3))
#  local samplename_num=$((($samplename_tmp +1)%2 + 1))
#  samplename=$(echo $samplenames | cut -d" " -f$samplename_num)
  local filter_num=$((($TASK_ID+2)/3))
  filter=$(echo $filters | cut -d" " -f$filter_num)
  extension=$(echo $filter_extensions | cut -d" " -f$filter_num)
}

get_filters $SGE_TASK_ID

samplenameS1=UC10164S1
sample_S1_vcf=${project_dir}/5Variants/1GATK_out_RefBias/${samplenameS1}/mapQ${mapQfilter}/${samplenameS1}${filter}.vcf
samplenameS2=UC10164S2
sample_S2_vcf=${project_dir}/5Variants/1GATK_out_RefBias/${samplenameS2}/mapQ${mapQfilter}/${samplenameS2}${filter}.vcf

outdir=${project_dir}/5Variants/2Merged_vcf_${extension}/mapQ${mapQfilter}
mkdir -p ${outdir}

#gzip and index all the vcfs
bcftools view --output-type z $sample_S1_vcf > ${TMPDIR}/sample_S1.vcf.gz
bcftools index ${TMPDIR}/sample_S1.vcf.gz
bcftools view --output-type z $sample_S2_vcf > ${TMPDIR}/sample_S2.vcf.gz
bcftools index ${TMPDIR}/sample_S2.vcf.gz
bcftools view --output-type z $outgroup_vcf > ${TMPDIR}/outgroup.vcf.gz
bcftools index ${TMPDIR}/outgroup.vcf.gz

bcftools merge ${TMPDIR}/sample_S1.vcf.gz ${TMPDIR}/sample_S2.vcf.gz ${TMPDIR}/outgroup.vcf.gz --output-type z > ${outdir}/outgroup_ancient.vcf.gz
bcftools index ${outdir}/outgroup_ancient.vcf.gz

