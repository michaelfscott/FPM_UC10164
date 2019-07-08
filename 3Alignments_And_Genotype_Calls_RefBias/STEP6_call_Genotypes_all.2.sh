#$ -cwd
#job name
#$ -N STEP6_call_Genotypes_all_2
#$ -o STEP6_call_Genotypes_all.2.log
#memory and runtime options 
#$ -l mem=30G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=70G

SGE_TASK_ID=8

#parameters
numthreads=1
java_memory_param=18
interval_padding_param=100
verbosity_param=ERROR

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
alternate_reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.alternate.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
species_altRef=Emmer_AltRef
species_altBam=Emmer_AltBam
read_type=Collapsed_Trimmed

#choose mapQ filter and samplename parameters across SGE_TASK_ID's
get_filters() {
  local TASK_ID=$1
  #in order to cycle through the chromosomes and subgenomes
  local mapQs="20 25 30"
  local samplenames="UC10164S1 UC10164S2"
  local filters=".bedfilter .bedfilter.noindels .bedfilter.posfilter.altref .bedfilter.noindels.posfilter.altBam .bedfilter.noindels.posfilter.altRef.altBam"
  local filter_extensions="bedfilter bedfilter_no_indels bedfilter_altref bedfilter_no_indels_altbam bedfilter_no_indels_altref_altbam"
  local mapQ_num=$((($TASK_ID+2)%3 + 1))
  mapQfilter=$(echo $mapQs | cut -d" " -f$mapQ_num)
  local samplename_tmp=$((($TASK_ID +2)/3))
  local samplename_num=$((($samplename_tmp +1)%2 + 1))
  samplename=$(echo $samplenames | cut -d" " -f$samplename_num)
  local filter_num=$((($TASK_ID+5)/6))
  filter=$(echo $filters | cut -d" " -f$filter_num)
  extension=$(echo $filter_extensions | cut -d" " -f$filter_num)
}

get_filters $SGE_TASK_ID

#will use merged bams with duplicates marked
input_bam=${project_dir}/4Alignments/3RefBias_Filters/${samplename}/mapQ${mapQfilter}/${samplename}${filter}.dedup.bam
#output dir
outdir=${project_dir}/5Variants/1GATK_out_RefBias/${samplename}
mkdir -p ${outdir}

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
--output ${outdir}/mapQ${mapQfilter}/${samplename}${filter}.vcf

