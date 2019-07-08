#$ -cwd
#$ -o STEP3_Count_Reads_Emmer_Trimmed.2.log
#$ -N STEP3_Count_Reads_Emmer_Trimmed.2
#memory and runtime options 
#$ -l mem=4G
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
read_type=Collapsed_Trimmed

samplenames="UC10164S1 UC10164S2"
for samplename in $samplenames; do 
echo $samplename

input_dir=${project_dir}/4Alignments/1Align_${species}_${read_type}/Merged/${samplename}
#get statistics for different filtering levels after dedup
input_bam=${input_dir}/${samplename}.dedup.bam
input_bam_full=${input_dir}/full_bam/${samplename}.bam

samtools bedcov ${known_variants}.bed \
${input_dir}/${samplename}.rl35.bam \
${input_dir}/${samplename}.dedup.rl35.bam \
${input_dir}/${samplename}.mapQ20.rl35.bam \
${input_dir}/${samplename}.dedup.mapQ20.rl35.bam \
${input_dir}/${samplename}.mapQ30.rl35.bam \
${input_dir}/${samplename}.dedup.mapQ30.rl35.bam > ${input_dir}/bedcov/SNP.rl35.bedcov

wait 

done


