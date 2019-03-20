#$ -cwd
#job name
#$ -N STEP3_Call_Genotypes
#$ -o STEP3_Call_Genotypes.log
#memory and runtime options 
#$ -l mem=30G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=40G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-2

#parameters
numthreads=1
java_memory_param=18
interval_padding_param=100
verbosity_param=ERROR

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
species=Emmer
read_type=Collapsed_Trimmed
#we will call variants at a set of known sites, specified in this vcf (split to match reference genome
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

#will use merged bams with duplicates marked
datadir=${project_dir}/4Alignments/1Align_${species}_${read_type}/Merged
#get samplename
ls -d ${datadir}/* > ${TMPDIR}/${species}_${read_type}_datadirs.txt
input_datadir=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/${species}_${read_type}_datadirs.txt)
samplename=$(echo $input_datadir | awk -F'/' '{print $(NF)}')
input_bam=${input_datadir}/${samplename}.dedup.bam
#output dir
outdir=${project_dir}/5Variants/1GATK_out/${samplename}
mkdir -p ${outdir}

for mapQfilter in 20 25 35; do 

#make a temporary dir for gatk to work in
mkdir -p ${TMPDIR}/${samplename}
mkdir -p ${outdir}/mapQ${mapQfilter}

echo "Haplotype Caller for ${samplename} mapQ ${mapQfilter} at $(date)"
#then perform HaplotypeCaller 
gatk HaplotypeCaller --java-options "-Xmx${java_memory_param}G" \
--reference ${reference} \
--TMP_DIR=${TMPDIR}/${samplename} \
--input $input_bam \
--alleles $known_variants \
--intervals ${known_variants} \
--interval-padding $interval_padding_param \
--minimum-mapping-quality ${mapQfilter} \
--read-filter ReadLengthReadFilter \
--min-read-length 35 \
--max-read-length 150 \
--verbosity $verbosity_param \
--genotyping-mode GENOTYPE_GIVEN_ALLELES \
--genotype-filtered-alleles true \
--output-mode EMIT_ALL_SITES \
--standard-min-confidence-threshold-for-calling 0 \
--output ${outdir}/mapQ${mapQfilter}/${samplename}_mapQ${mapQfilter}.vcf

done
