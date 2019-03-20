#$ -cwd
#job name
#$ -N STEP4_Merge_VCFs
#$ -o STEP4_Merge_VCFs.log
#memory and runtime options 
#$ -l mem=4G
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=20G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-3

#parameters
numthreads=1
max_depth_filter=30

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
#we will merge the variants with the set of modern emmer wheat variants
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

#repeat this across all of the 
mapQfilter_list="20 25 35"
mapQfilter=$(echo $mapQfilter_list | cut -d" " -f$SGE_TASK_ID) 

samplenameS1=UC10164S1
sample_S1_vcf=${project_dir}/5Variants/1GATK_out/${samplenameS1}/mapQ${mapQfilter}/${samplenameS1}_mapQ${mapQfilter}.vcf
samplenameS2=UC10164S2
sample_S2_vcf=${project_dir}/5Variants/1GATK_out/${samplenameS2}/mapQ${mapQfilter}/${samplenameS2}_mapQ${mapQfilter}.vcf

outdir=${project_dir}/5Variants/2Merged_vcf/mapQ${mapQfilter}
mkdir -p ${outdir}

#because the AD field clashes (R veersus .), we remove it so that we can merge the vcfs
bcftools annotate --remove FORMAT/AD --output-type z $sample_S1_vcf > ${TMPDIR}/sample_S1.vcf.gz
bcftools index ${TMPDIR}/sample_S1.vcf.gz
bcftools annotate --remove FORMAT/AD --output-type z $sample_S2_vcf > ${TMPDIR}/sample_S2.vcf.gz
bcftools index ${TMPDIR}/sample_S2.vcf.gz

bcftools merge ${TMPDIR}/sample_S1.vcf.gz ${TMPDIR}/sample_S2.vcf.gz ${known_variants}.gz --output-type z > ${outdir}/merged.vcf.gz

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


#we use a high confidence and low confidence set, with depth filters 2 and 4 on sample S1
for minDP_S1 in 2 4; do

mkdir -p ${outdir}/S1minDP${minDP_S1}

bcftools view ${outdir}/merged.vcf.gz \
--include "FORMAT/DP[0] >= ${minDP_S1}" \
--output-type v \
--output-file ${TMPDIR}/merged_S1minDP${minDP_S1}.vcf

vcftools --vcf ${TMPDIR}/merged_S1minDP${minDP_S1}.vcf \
--plink --chrom-map ${TMPDIR}/chromosome_map.txt \
--out ${TMPDIR}/merged_S1minDP${minDP_S1}

plink \
--file ${TMPDIR}/merged_S1minDP${minDP_S1} \
--recode A-transpose \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${outdir}/S1minDP${minDP_S1}/merged

done

#for comparison between samples, we will look at depth thresholds of 1 and 2
for minDP_S1 in 1 2; do
for minDP_S2 in 1 2; do

mkdir -p ${outdir}/S1minDP${minDP_S1}_S2minDP${minDP_S2}

bcftools view ${outdir}/merged.vcf.gz \
--include "FORMAT/DP[0] >= ${minDP_S1} && FORMAT/DP[1] >= ${minDP_S2}" \
--output-type v \
--output-file ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2}.vcf

vcftools --vcf ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2}.vcf \
--plink --chrom-map ${TMPDIR}/chromosome_map.txt \
--out ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2}

plink \
--file ${TMPDIR}/merged_S1minDP${minDP_S1}_S2minDP${minDP_S2} \
--recode A-transpose \
--reference-allele ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${outdir}/S1minDP${minDP_S1}_S2minDP${minDP_S2}/merged

done
done
