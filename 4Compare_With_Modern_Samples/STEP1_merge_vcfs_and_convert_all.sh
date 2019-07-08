#$ -cwd
#job name
#$ -N STEP1_merge_vcfs_and_convert_all
#$ -o STEP1_merge_vcfs_and_convert_all.log
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

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#we will merge the variants with the set of modern emmer wheat variants
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

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

outdir=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}
mkdir -p ${outdir}
#here we read in the ancient and outgroup genotypes from vcf. If the genotyping step is skipped, the file position should be specified here
ancient_outgroup_variants=${project_dir}/5Variants/2Merged_vcf_${extension}/mapQ${mapQfilter}/outgroup_ancient.vcf.gz

#because the AD field clashes (R veersus .), we remove it so that we can merge the vcfs
bcftools annotate --remove FORMAT/AD --output-type z $ancient_outgroup_variants > ${TMPDIR}/ancient_outgroup.vcf.gz
bcftools index ${TMPDIR}/ancient_outgroup.vcf.gz

#merge with modern variants
bcftools merge ${TMPDIR}/ancient_outgroup.vcf.gz ${known_variants}.gz --output-type z > ${outdir}/merged.vcf.gz

#for plink, we speficy the allele to count so that we count alternate alleles:
awk -v OFS='\t' '{print $1":"$2,$5}' ${known_variants} | \
grep -v "^#" > ${TMPDIR}/reference_alleles.txt
#should also specify the chromosome map so that the chromosomes are numbered sequentially:
awk -v OFS='\t' '{print $1, NR}' ${reference}.fai > ${TMPDIR}/chromosome_map.txt

#output the full traw plink file
vcftools --gzvcf ${outdir}/merged.vcf.gz \
--plink --chrom-map ${TMPDIR}/chromosome_map.txt \
--out ${TMPDIR}/merged
plink \
--file ${TMPDIR}/merged \
--recode A-transpose \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${outdir}/merged
#also output only selected regions
plink --noweb \
--recode A-transpose \
--extract range selected_regions.txt \
--chr-set 30 \
--file ${TMPDIR}/merged \
--out ${outdir}/merged.regions

#we also wish to filter on depth and then output plink files for analysis
#we use a high confidence and low confidence set, with depth filters 2 and 4 on sample S1
for minDP_S1 in 2 4; do

mkdir -p ${outdir}/S1minDP${minDP_S1}

#filter S1 depth
bcftools view ${outdir}/merged.vcf.gz \
--include "FORMAT/DP[0] >= ${minDP_S1} && FORMAT/DP[0] < ${max_depth_filter}" \
--output-type v \
--output-file ${TMPDIR}/merged_S1minDP${minDP_S1}.vcf
#convert for plink
vcftools --vcf ${TMPDIR}/merged_S1minDP${minDP_S1}.vcf \
--plink --chrom-map ${TMPDIR}/chromosome_map.txt \
--out ${TMPDIR}/merged_S1minDP${minDP_S1}
#output traw
plink \
--file ${TMPDIR}/merged_S1minDP${minDP_S1} \
--recode A-transpose \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${outdir}/S1minDP${minDP_S1}/merged
#output bfiles 
plink \
--file ${TMPDIR}/merged_S1minDP${minDP_S1} \
--make-bed \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${outdir}/S1minDP${minDP_S1}/merged

done

#for comparison between samples, we will look at depth thresholds of 1 and 2
for minDP_S1 in 1 2; do
for minDP_S2 in 1 2; do

mkdir -p ${outdir}/S1minDP${minDP_S1}_S2minDP${minDP_S2}
#filter on depth
bcftools view ${outdir}/merged.vcf.gz \
--include "FORMAT/DP[0] >= ${minDP_S1} && FORMAT/DP[1] >= ${minDP_S2} && FORMAT/DP[0] < ${max_depth_filter} && FORMAT/DP[1] < ${max_depth_filter}" \
--output-type v \
--output-file ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2}.vcf
#convert for plink
vcftools --vcf ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2}.vcf \
--plink --chrom-map ${TMPDIR}/chromosome_map.txt \
--out ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2}
#output traw
plink \
--file ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2} \
--recode A-transpose \
--reference-allele ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${outdir}/S1minDP${minDP_S1}_S2minDP${minDP_S2}/merged

done
done
