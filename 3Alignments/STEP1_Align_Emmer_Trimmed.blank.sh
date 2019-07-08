#$ -cwd
#job name
#$ -N STEP1_Align_Emmer_blank
#$ -o STEP1_Align_Emmer.blank.log
#memory and runtime options 
#$ -l mem=8G
#Number of threads per job, memory given above is per-thread
#$ -pe smp 4
#$ -l h_rt=24:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=60G
#request number of threads (each has different SGE_TASK_ID variable)

#export variables
set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
species=Emmer
read_type=Collapsed_Trimmed
datadir=${project_dir}/3AdapterRemoval/3Adapters_Removed_And_Trimmed

#parameters
numthreads=4
bwa_seedlength=16500
bwa_editdist=0.01
bwa_gapopens=2

#we look in the datadir for the single end reads after adapter removal
#for the library blank, we can combine them because there are very few reads
fastq=${TMPDIR}/Sel_II_S2.fastq
cat ${datadir}/blank/Sel_II_S2/*.collapsed.not_truncated.fastq.gz > $fastq

samplename=blank
prefix=${samplename}
experimentname=Sel_II_S2

#create output directory
outdir=${project_dir}/4Alignments/1Align_${species}_${read_type}/Merged/${samplename}/
mkdir -p ${outdir}

#make a temporary directory for samtools sort
mkdir -p ${TMPDIR}/sort_${prefix}

echo "$SGE_TASK_ID aligning SE reads $prefix at $(date)"
bwa aln -t $numthreads -l $bwa_seedlength -n $bwa_editdist -o $bwa_gapopens $reference ${fastq} > ${TMPDIR}/${prefix}_sa.sai 
echo "$SGE_TASK_ID samse and samtools sort for $prefix at $(date)"
bwa samse $reference ${TMPDIR}/${prefix}_sa.sai ${fastq} | \
samtools sort -@ $numthreads -T ${TMPDIR}/sort_${prefix}/${prefix} -o ${TMPDIR}/${prefix}.sorted.bam -
echo "$SGE_TASK_ID indexing ${prefix} at $(date)"
samtools index ${TMPDIR}/${prefix}.sorted.bam

echo "$SGE_TASK_ID copying at $(date)"
cp ${TMPDIR}/${prefix}.sorted.bam* --target-directory=${outdir}/

samtools stats ${TMPDIR}/${}.bam > ${outdir}/stats/${samplename}.stats
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

