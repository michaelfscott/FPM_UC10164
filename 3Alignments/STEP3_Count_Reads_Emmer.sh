#$ -cwd
#$ -o STEP3_Count_Reads_Emmer.log
#$ -N STEP3_Count_Reads_Emmer
#memory and runtime options 
#$ -l mem=3G
#$ -pe smp 6
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#request TMPDIR space
#$ -l tmpfs=10G
#request number of threads (each has different SGE_TASK_ID variable)

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

#get statistics for different filtering levels

#the following specifies whether we use the trimmed or non-trimmed reads and which reference genome to use
reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
species=Emmer
read_type=Collapsed

samplenames="UC10164S1 UC10164S2"
for samplename in $samplenames; do 
echo $samplename

input_dir=${project_dir}/4Alignments/1Align_${species}_${read_type}/Merged/${samplename}
#get statistics for different filtering levels after dedup
input_bam=${input_dir}/${samplename}.dedup.bam
input_bam_full=${input_dir}/full_bam/${samplename}.bam

samtools view -h $input_bam | awk 'length($10) > 34 || $1 ~ /^@/' | samtools view -b > ${input_dir}/${samplename}.dedup.rl35.bam &
samtools view -h $input_bam_full | awk 'length($10) > 34 || $1 ~ /^@/' | samtools view -b > ${input_dir}/${samplename}.rl35.bam &

for mapQfilter in 20 30; do
samtools view -q $mapQfilter -b $input_bam > ${input_dir}/${samplename}.dedup.mapQ${mapQfilter}.bam &
samtools view -q $mapQfilter -b $input_bam_full > ${input_dir}/${samplename}.mapQ${mapQfilter}.bam &
samtools view -q $mapQfilter -h $input_bam | awk 'length($10) > 34 || $1 ~ /^@/' | samtools view -b > ${input_dir}/${samplename}.dedup.mapQ${mapQfilter}.rl35.bam &
samtools view -q $mapQfilter -h $input_bam_full | awk 'length($10) > 34 || $1 ~ /^@/' | samtools view -b > ${input_dir}/${samplename}.mapQ${mapQfilter}.rl35.bam &
done 
wait

samtools index ${input_dir}/${samplename}.dedup.rl35.bam &
samtools index ${input_dir}/${samplename}.rl35.bam &
for mapQfilter in 20 30; do
samtools index ${input_dir}/${samplename}.dedup.mapQ${mapQfilter}.bam &
samtools index ${input_dir}/${samplename}.mapQ${mapQfilter}.bam &
samtools index ${input_dir}/${samplename}.dedup.mapQ${mapQfilter}.rl35.bam &
samtools index ${input_dir}/${samplename}.mapQ${mapQfilter}.rl35.bam &
done
wait

samtools stats ${input_dir}/${samplename}.dedup.rl35.bam > ${input_dir}/stats/${samplename}.dedup.rl35.stats &
samtools stats ${input_dir}/${samplename}.rl35.bam > ${input_dir}/stats/${samplename}.rl35.stats &

for mapQfilter in 20 30; do
samtools stats ${input_dir}/${samplename}.dedup.mapQ${mapQfilter}.bam > ${input_dir}/stats/${samplename}.dedup.mapQ${mapQfilter}.stats &
samtools stats ${input_dir}/${samplename}.mapQ${mapQfilter}.bam > ${input_dir}/stats/${samplename}.mapQ${mapQfilter}.stats &
samtools stats ${input_dir}/${samplename}.dedup.mapQ${mapQfilter}.rl35.bam > ${input_dir}/stats/${samplename}.dedup.mapQ${mapQfilter}.rl35.stats &
samtools stats ${input_dir}/${samplename}.mapQ${mapQfilter}.rl35.bam > ${input_dir}/stats/${samplename}.mapQ${mapQfilter}.rl35.stats &
done

awk -v OFS='\t' '{print $1, 0, $2-1}' ${reference}.fai > ${TMPDIR}/all_chr.bed
samtools bedcov ${TMPDIR}/all_chr.bed \
${input_dir}/${samplename}.rl35.bam \
${input_dir}/${samplename}.dedup.rl35.bam \
${input_dir}/${samplename}.mapQ20.rl35.bam \
${input_dir}/${samplename}.dedup.mapQ20.rl35.bam \
${input_dir}/${samplename}.mapQ30.rl35.bam \
${input_dir}/${samplename}.dedup.mapQ30.rl35.bam > ${input_dir}/bedcov/chromosomal.rl35.bedcov

samtools bedcov ${known_variants}.bed \
${input_dir}/${samplename}.rl35.bam \
${input_dir}/${samplename}.dedup.rl35.bam \
${input_dir}/${samplename}.mapQ20.rl35.bam \
${input_dir}/${samplename}.dedup.mapQ20.rl35.bam \
${input_dir}/${samplename}.mapQ30.rl35.bam \
${input_dir}/${samplename}.dedup.mapQ30.rl35.bam > ${input_dir}/bedcov/SNP.rl35.bedcov

wait 

done


