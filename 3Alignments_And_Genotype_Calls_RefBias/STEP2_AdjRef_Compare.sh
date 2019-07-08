#$ -cwd
#job name
#$ -N STEP2_AdjRef_Compare
#$ -o STEP2_AdjRef_Compare.log
#memory and runtime options 
#$ -l mem=5G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=80G
#request number of threads (each has different SGE_TASK_ID variable)
#$ -t 1-89

#parameters
numthreads=1

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh
module load python2/recommended

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
alternate_reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.alternate.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
species_altRef=Emmer_AltRef
read_type=Collapsed_Trimmed

#get samplename and fastq file to work with
datadir=${project_dir}/4Alignments/2Align_${species}_${read_type}/SplitReads
ls ${datadir}/*/*.sorted.bam > ${TMPDIR}/${species}_${read_type}_bams.txt
datadir_alt=${project_dir}/4Alignments/2Align_${species_altRef}_${read_type}/SplitReads
ls ${datadir_alt}/*/*.sorted.bam > ${TMPDIR}/${species_altRef}_${read_type}_bams.txt

original=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/${species}_${read_type}_bams.txt )
alternate=$(sed "${SGE_TASK_ID}q;d" ${TMPDIR}/${species_altRef}_${read_type}_bams.txt )

prefix=$( basename $original .trimmed.sorted.bam)
samplename=$( echo $original | awk -v FS="/" '{print $(NF-1)}' )

#make bed file from vcf
awk -v FS='\t' -v OFS='\t' '/^[^#]/ { $2 = $2 - 1; print $1, $2, $2 }' $known_variants > ${TMPDIR}/split_SNPs.bed
bedfile=${TMPDIR}/split_SNPs.bed

for mapQfilter in 20 25 30; do 
#output_directory
bedfilter_dir=${project_dir}/4Alignments/2Align_${species_altRef}_${read_type}/bedfilter/${samplename}/mapQ${mapQfilter}
mkdir -p $bedfilter_dir

#calmd to get full comparison with the reference. 
#filters to regions of interest and uses a quality threshold of 30
samtools calmd ${original} ${reference} 2> ${TMPDIR}/orig_calmd.log | samtools view -L $bedfile -h -b -q $mapQfilter - > ${TMPDIR}/bedfilter.bam 
samtools index ${TMPDIR}/bedfilter.bam
samtools calmd ${alternate} ${alternate_reference} 2> ${TMPDIR}/alt_calmd.log | samtools view -L $bedfile -h -b -q $mapQfilter - > ${TMPDIR}/bedfilter.altref.bam 
samtools index ${TMPDIR}/bedfilter.altref.bam

#posfilter for bedfilter.bam, leaves alignment against original reference genome in the bam file
samtools view -@ $numthreads -h ${TMPDIR}/bedfilter.bam | \
python2 ${REFBIAS_PATH}/filter_sam_startpos_dict.mapQ${mapQfilter}.py ${TMPDIR}/bedfilter.altref.bam | \
samtools view -@ $numthreads -bh - > ${TMPDIR}/bedfilter.posfilter.altref.bam
samtools index ${TMPDIR}/bedfilter.posfilter.altref.bam

for fb in bedfilter bedfilter.altref bedfilter.posfilter.altref 
do
mv ${TMPDIR}/$fb.bam ${bedfilter_dir}/$prefix.$fb.bam
mv ${TMPDIR}/$fb.bam.bai ${bedfilter_dir}/$prefix.$fb.bam.bai
done

done
