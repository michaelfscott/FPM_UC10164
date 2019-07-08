#$ -cwd
#job name
#$ -N STEP4_Align_Emmer_Trimmed_AdjBam
#$ -o STEP4_Align_Emmer_Trimmed_AdjBam.log
#memory and runtime options 
#$ -l mem=8G
#$ -pe smp 4
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=80G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-6

#parameters
numthreads=4
bwa_seedlength=16500
bwa_editdist=0.01
bwa_gapopens=2
java_memory_param=18

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
#need to create a SNPlist 
awk -v FS='\t' -v OFS='\t' '/^[^#]/ {print $1, $2, $4, $5 }' $known_variants > ${TMPDIR}/snplist.txt
#make bed file from vcf
bedfile=${TMPDIR}/split_SNPs.bed
awk -v FS='\t' -v OFS='\t' '/^[^#]/ { $2 = $2 - 1; print $1, $2, $2 }' $known_variants > $bedfile
#also need chromosome names because they are non-numeric:
awk '{print $1}' $reference.fai > ${TMPDIR}/chromosome_names.txt

##############
#bedfilter file input
bedfilter_Merged_dir=${project_dir}/4Alignments/2Align_${species_altRef}_${read_type}/Merged_bedfilter/${samplename}/mapQ${mapQfilter}
cp ${bedfilter_Merged_dir}/${samplename}.bedfilter.dedup.bam ${TMPDIR}/bedfilter.bam
#output dir
bedfilter_Merged_altBam_dir=${project_dir}/4Alignments/2Align_${species_altBam}_${read_type}/Merged_bedfilter/${samplename}/mapQ${mapQfilter}
refbias_filter_dir=${project_dir}/4Alignments/3RefBias_Filters/${samplename}/mapQ${mapQfilter}
mkdir -p ${bedfilter_Merged_altBam_dir}
mkdir -p ${refbias_filter_dir}

#the python script that modifies reads does not allow indels
#therefore we also remove indels from the original bam for comparison by using the custom python script below
samtools calmd ${TMPDIR}/bedfilter.bam $reference | \
python2 ${REFBIAS_PATH}/modify_read_filter_indels.py ${TMPDIR}/snplist.txt ${TMPDIR}/chromosome_names.txt | \
samtools view -@ $numthreads -b -h - > ${TMPDIR}/bedfilter.noindels.bam
samtools index ${TMPDIR}/bedfilter.noindels.bam

#adjust reads from bam
samtools calmd ${TMPDIR}/bedfilter.bam $reference | \
python2 ${REFBIAS_PATH}/modify_read_alternative.py ${TMPDIR}/snplist.txt ${TMPDIR}/chromosome_names.txt | \
samtools view -@ $numthreads -b -h - > ${TMPDIR}/modified.noindels.bam
samtools index ${TMPDIR}/bedfilter.noindels.bam

#re-align
bwa aln -t $numthreads -l $bwa_seedlength -n $bwa_editdist -o $bwa_gapopens $reference -b ${TMPDIR}/modified.noindels.bam | \
bwa samse ${reference} - ${TMPDIR}/modified.noindels.bam  | \
samtools view -@ $numthreads -F 4 -h -b - | \
samtools sort -o ${TMPDIR}/modified.noindels.remap.bam -
samtools index ${TMPDIR}/modified.noindels.remap.bam

#get MD 
samtools calmd ${TMPDIR}/modified.noindels.remap.bam $reference | samtools view -@ $numthreads -bh - > ${TMPDIR}/modified.noindels.remap.MD.bam
samtools index ${TMPDIR}/modified.noindels.remap.MD.bam

#filter based on alignment position overlap
samtools view -@ $numthreads -h ${TMPDIR}/bedfilter.noindels.bam | \
python2 ${REFBIAS_PATH}/filter_sam_startpos_dict.mapQ${mapQfilter}.py ${TMPDIR}/modified.noindels.remap.MD.bam | \
samtools view -@ $numthreads -bh - > ${TMPDIR}/bedfilter.noindels.posfilter.altBam.bam
samtools index ${TMPDIR}/bedfilter.noindels.posfilter.altBam.bam

#output results
for fb in bedfilter.noindels modified.noindels modified.noindels.remap bedfilter.noindels.posfilter.altBam
do
cp ${TMPDIR}/$fb.bam ${bedfilter_Merged_altBam_dir}/${samplename}.$fb.bam
cp ${TMPDIR}/$fb.bam.bai ${bedfilter_Merged_altBam_dir}/${samplename}.$fb.bam.bai
done

### dedup 
#output dedup metrics to file in this dir
mkdir -p ${bedfilter_Merged_altBam_dir}/dedup_metrics/

#cycle through the filters
for filter in bedfilter.noindels modified.noindels modified.noindels.remap bedfilter.noindels.posfilter.altBam; do

#make a temporary directory for dedup
mkdir -p ${TMPDIR}/${samplename}${filter}
#mark Duplicates and index
echo "$SGE_TASK_ID mark Duplicates for ${samplename} ${species} ${read_type} at $(date)"
gatk MarkDuplicates --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=${TMPDIR}/${samplename}${filter} \
--INPUT ${TMPDIR}/${filter}.bam \
--METRICS_FILE ${bedfilter_Merged_altBam_dir}/dedup_metrics/${samplename}${filter}.mapQ${mapQfilter}.dedup_metrics.txt \
--OUTPUT ${TMPDIR}/${filter}.dedup.bam
echo "$SGE_TASK_ID indexing merged dedup bam ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/${filter}.dedup.bam

cp ${TMPDIR}/${filter}.dedup.bam ${bedfilter_Merged_altBam_dir}/${samplename}.${filter}.dedup.bam
cp ${TMPDIR}/${filter}.dedup.bam.bai ${bedfilter_Merged_altBam_dir}/${samplename}.${filter}.dedup.bam.bai

done

for filter in bedfilter.noindels bedfilter.noindels.posfilter.altBam; do
cp ${TMPDIR}/${filter}.dedup.bam ${refbias_filter_dir}/${samplename}.${filter}.dedup.bam
cp ${TMPDIR}/${filter}.dedup.bam.bai ${refbias_filter_dir}/${samplename}.${filter}.dedup.bam.bai
done


