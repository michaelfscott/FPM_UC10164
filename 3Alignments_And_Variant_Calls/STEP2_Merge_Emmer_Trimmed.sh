#$ -cwd
#job name
#$ -N STEP2_Merge_Emmer_Trimmed
#$ -o STEP2_Merge_Emmer_Trimmed.log
#memory and runtime options 
#$ -l mem=30G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=80G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-2

#parameters
numthreads=1
java_memory_param=18

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
read_type=Collapsed_Trimmed

datadir=${project_dir}/4Alignments/1Align_${species}_${read_type}/SplitReads
ls -d ${datadir}/* > ${TMPDIR}/${species}_${read_type}_datadirs.txt
input_datadir=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/${species}_${read_type}_datadirs.txt)
samplename=$(echo $input_datadir | awk -F'/' '{print $(NF)}')
outdir=${project_dir}/4Alignments/1Align_${species}_${read_type}/Merged/${samplename}
mkdir -p ${outdir}

#make a list of bams to merge
bamlist=${TMPDIR}/${species}_${read_type}_${samplename}_bams.txt
ls ${input_datadir}/*.sorted.bam > $bamlist

#the c and p options are so that RG and PG tags with colliding IDs are merged rather than adding a suffix to differentiate them. 
echo "$SGE_TASK_ID merging ${samplename} ${species} ${read_type} at $(date)"
samtools merge -c -p -b $bamlist ${TMPDIR}/${samplename}.bam 
echo "$SGE_TASK_ID indexing ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/${samplename}.bam

#get alignment statistics before marking duplicates (and filtering unmapped reads)
mkdir -p ${outdir}/stats
echo "$SGE_TASK_ID alignment statistics ${samplename} ${species} ${read_type} at $(date)"
samtools stats ${TMPDIR}/${samplename}.bam > ${outdir}/stats/${samplename}.stats
grep "^RL" ${outdir}/stats/${samplename}.stats | cut -f 2- > ${outdir}/stats/read_lengths.txt

#for picard markDuplicates, need to remove unmapped reads first. 
echo "filtering ${samplename} ${species} ${read_type} to mapped reads only at $(date)"
samtools view -b -h -F 4 ${TMPDIR}/${samplename}.bam > ${TMPDIR}/${samplename}.mapped.bam
samtools index ${TMPDIR}/${samplename}.mapped.bam

#make a temporary directory for dedup
mkdir -p ${TMPDIR}/${samplename}
#output dedup metrics to file in this dir
mkdir -p ${outdir}/dedup_metrics/

#mark Duplicates and index
echo "$SGE_TASK_ID mark Duplicates for ${samplename} ${species} ${read_type} at $(date)"
gatk MarkDuplicates --java-options "-Xmx${java_memory_param}G" \
--TMP_DIR=${TMPDIR}/${samplename} \
--INPUT ${TMPDIR}/${samplename}.mapped.bam \
--METRICS_FILE ${outdir}/dedup_metrics/${samplename}.dedup_metrics.txt \
--OUTPUT ${TMPDIR}/${samplename}.dedup.bam 
echo "$SGE_TASK_ID indexing merged dedup bam ${samplename} ${species} ${read_type} at $(date)"
samtools index ${TMPDIR}/${samplename}.dedup.bam 

#copy out dedup file
echo "$SGE_TASK_ID copying at $(date)"
cp ${TMPDIR}/${samplename}.dedup.bam* --target-directory=${outdir}/

#calculate statistics and mapDamage profile
#first make output folders for results
mkdir -p ${outdir}/stats
mkdir -p ${outdir}/map_damage
echo "$SGE_TASK_ID alignment statistics ${samplename} ${species} ${read_type} at $(date)"
samtools stats ${TMPDIR}/${samplename}.dedup.bam > ${outdir}/stats/${samplename}.dedup.stats 
echo "$SGE_TASK_ID mapdamage for ${samplename} ${species} ${read_type} at $(date)"
mapDamage -i ${TMPDIR}/${samplename}.dedup.bam -r $reference --folder ${outdir}/map_damage

#calculate coverage over chromosomes
mkdir -p ${outdir}/bedcov
#first get the start and end positions of each chromosome in a bed file (0 based)
awk -v OFS='\t' '{print $1, 0, $2-1}' ${reference}.fai > ${TMPDIR}/all_chr.bed
for mapQfilter in 0 1 10 20 25 30 35; do
samtools bedcov -Q $mapQfilter ${TMPDIR}/all_chr.bed ${TMPDIR}/${samplename}.bam ${TMPDIR}/${samplename}.dedup.bam > ${outdir}/bedcov/chromosomal.mapQ${mapQfilter}.bedcov
done
#for Emmer wheat, we will also specifically check coverage of the regions where SNPs will be called
for mapQfilter in 0 1 10 20 25 30 35; do
samtools bedcov -Q $mapQfilter ${known_variants}.bed ${TMPDIR}/${samplename}.bam ${TMPDIR}/${samplename}.dedup.bam > ${outdir}/bedcov/SNP.mapQ${mapQfilter}.bedcov
done


