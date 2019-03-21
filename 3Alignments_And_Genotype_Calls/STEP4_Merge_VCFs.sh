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

#repeat this across all of the mapQ's
mapQfilter_list="20 25 35"
mapQfilter=$(echo $mapQfilter_list | cut -d" " -f$SGE_TASK_ID) 

samplenameS1=UC10164S1
sample_S1_vcf=${project_dir}/5Variants/1GATK_out/${samplenameS1}/mapQ${mapQfilter}/${samplenameS1}_mapQ${mapQfilter}.vcf
samplenameS2=UC10164S2
sample_S2_vcf=${project_dir}/5Variants/1GATK_out/${samplenameS2}/mapQ${mapQfilter}/${samplenameS2}_mapQ${mapQfilter}.vcf

outdir=${project_dir}/5Variants/2Merged_vcf/mapQ${mapQfilter}
mkdir -p ${outdir}

#gzip and index all the vcfs
bcftools view --output-type z $sample_S1_vcf > ${TMPDIR}/sample_S1.vcf.gz
bcftools index ${TMPDIR}/sample_S1.vcf.gz
bcftools view --output-type z $sample_S2_vcf > ${TMPDIR}/sample_S2.vcf.gz
bcftools index ${TMPDIR}/sample_S2.vcf.gz
bcftools view --output-type z $outgroup_vcf > ${TMPDIR}/outgroup.vcf.gz
bcftools index ${TMPDIR}/outgroup.vcf.gz

bcftools merge ${TMPDIR}/sample_S1.vcf.gz ${TMPDIR}/sample_S2.vcf.gz ${TMPDIR}/outgroup.vcf.gz --output-type z > ${outdir}/outgroup_ancient.vcf.gz



