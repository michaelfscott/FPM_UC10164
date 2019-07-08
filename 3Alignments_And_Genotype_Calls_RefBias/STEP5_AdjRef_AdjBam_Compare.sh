#$ -cwd
#job name
#$ -N STEP5_AdjRef_AdjBam_Compare
#$ -o STEP5_AdjRef_AdjBam_Compare.log
#memory and runtime options 
#$ -l mem=5G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=80G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-6

#parameters
numthreads=1
java_memory_param=2

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh
module load python2/recommended

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
alternate_reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.alternate.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
species_altRef=Emmer_AltRef
species_altBam=Emmer_AltBam
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

##############

altRef_posfilter=${project_dir}/4Alignments/2Align_${species_altRef}_${read_type}/Merged_bedfilter/${samplename}/mapQ${mapQfilter}/${samplename}.bedfilter.posfilter.altref.dedup.bam
altBam_posfilter=${project_dir}/4Alignments/2Align_${species_altBam}_${read_type}/Merged_bedfilter/${samplename}/mapQ${mapQfilter}/${samplename}.bedfilter.noindels.posfilter.altBam.dedup.bam

#output dir
refbias_filter_dir=${project_dir}/4Alignments/3RefBias_Filters/${samplename}/mapQ${mapQfilter}
mkdir -p ${refbias_filter_dir}

#calmd to get full comparison with the reference. 
#filters to regions of interest and uses a quality threshold of mapQfilter
samtools calmd ${altRef_posfilter} ${reference} 2> ${TMPDIR}/altRef_calmd.log | samtools view -h -b -q $mapQfilter - > ${TMPDIR}/bedfilter.posfilter.altref.dedup.MD.bam
samtools index ${TMPDIR}/bedfilter.posfilter.altref.dedup.MD.bam
samtools calmd ${altBam_posfilter} ${reference} 2> ${TMPDIR}/altBam_calmd.log | samtools view -h -b -q $mapQfilter - > ${TMPDIR}/bedfilter.noindels.posfilter.altBam.dedup.MD.bam
samtools index ${TMPDIR}/bedfilter.noindels.posfilter.altBam.dedup.MD.bam

#posfilter for bedfilter.bam, leaves alignment against original reference genome in the bam file
samtools view -@ $numthreads -h ${TMPDIR}/bedfilter.posfilter.altref.dedup.MD.bam | \
python2 ${REFBIAS_PATH}/filter_sam_startpos_dict.mapQ${mapQfilter}.py ${TMPDIR}/bedfilter.noindels.posfilter.altBam.dedup.MD.bam | \
samtools view -@ $numthreads -bh - > ${TMPDIR}/bedfilter.noindels.posfilter.altRef.altBam.bam
samtools index ${TMPDIR}/bedfilter.noindels.posfilter.altRef.altBam.bam

#make a temporary directory for dedup
mkdir -p ${TMPDIR}/${samplename}
mkdir -p ${refbias_filter_dir}/dedup_metrics
#mark Duplicates and index
echo "$SGE_TASK_ID mark Duplicates for ${samplename} ${species} ${read_type} at $(date)"
gatk MarkDuplicates --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=${TMPDIR}/${samplename} \
--INPUT ${TMPDIR}/bedfilter.noindels.posfilter.altRef.altBam.bam \
--METRICS_FILE ${refbias_filter_dir}/dedup_metrics/${samplename}.bedfilter.noindels.posfilter.altRef.altBam.dedup_metrics.mapQ${mapQfilter}.dedup_metrics.txt \
--OUTPUT ${TMPDIR}/bedfilter.noindels.posfilter.altRef.altBam.dedup.bam
echo "$SGE_TASK_ID indexing merged dedup bam ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/bedfilter.noindels.posfilter.altRef.altBam.dedup.bam


fb=bedfilter.noindels.posfilter.altRef.altBam.dedup
mv ${TMPDIR}/$fb.bam ${refbias_filter_dir}/${samplename}.$fb.bam
mv ${TMPDIR}/$fb.bam.bai ${refbias_filter_dir}/${samplename}.$fb.bam.bai

