#$ -cwd
#job name
#$ -N STEP3_Merge_Emmer_Trimmed_AdjRef
#$ -o STEP3_Merge_Emmer_Trimmed_AdjRef.log
#memory and runtime options 
#$ -l mem=10G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=20G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-6

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#parameters
numthreads=1
java_memory_param=5

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
alternate_reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.alternate.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
species_altRef=Emmer_AltRef
read_type=Collapsed_Trimmed

get_samplename_mapQfilter() {
  local TASK_ID=$1
  #in order to cycle through the chromosomes and subgenomes
  local mapQs="20 25 30"
#  local samplenames="UC10164S1 UC10164S2"
  local samplenames="UC10164S1 UC10164S2"
  local mapQ_num=$((($TASK_ID+2)%3 + 1))
  mapQfilter=$(echo $mapQs | cut -d" " -f$mapQ_num)
#  local samplename_tmp=$((($TASK_ID +2)/3))
#  local samplename_num=$((($samplename_tmp +1)%2 + 1))
#  samplename=$(echo $samplenames | cut -d" " -f$samplename_num)
  local samplename_num=$((($TASK_ID+2)/3))
  samplename=$(echo $samplenames | cut -d" " -f$samplename_num)
}

get_samplename_mapQfilter $SGE_TASK_ID

input_datadir=${project_dir}/4Alignments/2Align_${species_altRef}_${read_type}/bedfilter/${samplename}/mapQ${mapQfilter}
bedfilter_Merged_dir=${project_dir}/4Alignments/2Align_${species_altRef}_${read_type}/Merged_bedfilter/${samplename}/mapQ${mapQfilter}
refbias_filter_dir=${project_dir}/4Alignments/3RefBias_Filters/${samplename}/mapQ${mapQfilter}
mkdir -p ${bedfilter_Merged_dir}
mkdir -p ${refbias_filter_dir}

for filter in .bedfilter .bedfilter.posfilter.altref; do

#make a list of bams to merge
bamlist=${TMPDIR}/${species}_${read_type}_${samplename}bams${filter}.txt
ls ${input_datadir}/*${filter}.bam > $bamlist

#the c and p options are so that RG and PG tags with colliding IDs are merged rather than adding a suffix to differentiate them. 
echo "$SGE_TASK_ID merging ${samplename} ${species} ${read_type} at $(date)"
samtools merge -c -p -b $bamlist ${TMPDIR}/${samplename}${filter}.bam 
echo "$SGE_TASK_ID indexing ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/${samplename}${filter}.bam

#for picard markDuplicates, need to remove unmapped reads first. 
echo "filtering ${samplename} ${species} ${read_type} to mapped reads only at $(date)"
samtools view -b -h -F 4 ${TMPDIR}/${samplename}${filter}.bam > ${TMPDIR}/${samplename}${filter}.mapped.bam
samtools index ${TMPDIR}/${samplename}${filter}.mapped.bam

#make a temporary directory for dedup
mkdir -p ${TMPDIR}/${samplename}${filter}
#output dedup metrics to file in this dir
mkdir -p ${bedfilter_Merged_dir}/dedup_metrics/

#mark Duplicates and index
echo "$SGE_TASK_ID mark Duplicates for ${samplename} ${species} ${read_type} at $(date)"
gatk MarkDuplicates --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=${TMPDIR}/${samplename}${filter} \
--INPUT ${TMPDIR}/${samplename}${filter}.mapped.bam \
--METRICS_FILE ${bedfilter_Merged_dir}/dedup_metrics/${samplename}${filter}.mapQ${mapQfilter}.dedup_metrics.txt \
--OUTPUT ${TMPDIR}/${samplename}${filter}.dedup.bam 
echo "$SGE_TASK_ID indexing merged dedup bam ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/${samplename}${filter}.dedup.bam 

#copy out dedup file
echo "$SGE_TASK_ID copying at $(date)"
cp ${TMPDIR}/${samplename}${filter}.bam --target-directory=${bedfilter_Merged_dir}/
cp ${TMPDIR}/${samplename}${filter}.dedup.bam* --target-directory=${bedfilter_Merged_dir}/
cp ${TMPDIR}/${samplename}${filter}.dedup.bam* --target-directory=${refbias_filter_dir}/

done

