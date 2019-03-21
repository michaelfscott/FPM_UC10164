#job name
#$ -N STEP1_Split_Emmer_Reference
#$ -o STEP1_Split_Emmer_Reference.log
# Set the working directory 
#$ -cwd
#memory and runtime options 
#$ -l mem=40G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=30G

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

java_memory_param=10

input_reference_genome=${emmer_wheat_reference_original}
output_dir=${project_dir}/2Reference_Genomes/Zavitan_v2_split
output_file=${output_dir}/151210_zavitan_v2_pseudomolecules_split.fasta
split_fasta_dir=${TMPDIR}/separate_chromosomes

mkdir -p $output_dir
mkdir -p $split_fasta_dir

#index the input reference genome
echo "samtools faidx emmer reference at $(date)"
samtools faidx $input_reference_genome

#The fai file has the names and lengths of all the chromosomes
awk '{print $1 "\t" $2}' ${input_reference_genome}.fai > ${TMPDIR}/chr_names_lengths.txt

#for each chromosome, split into two halves. ${chr}_1 has the positions from 1-400mb, ${chr}_2 has the remaining positions.
while read line ; do
  chr=$(echo $line | awk '{print $1}')
  length=$(echo $line | awk '{print $2}')
  echo $chr
  echo $length
  samtools faidx $input_reference_genome ${chr}:1-400000000 > ${split_fasta_dir}/${chr}_1.fasta 
  samtools faidx $input_reference_genome ${chr}:400000001-${length} > ${split_fasta_dir}/${chr}_2.fasta 
done < ${TMPDIR}/chr_names_lengths.txt

#join all split chromosome halves into one fasta file to get new split reference genome
cat ${split_fasta_dir}/*.fasta > ${output_file}

#index reference genome
echo "bwa indexing reference at $(date)"
bwa index ${output_file}
echo "samtools faidx reference at $(date)"
samtools faidx ${output_file}
echo "GATK indexing reference at $(date)"
gatk CreateSequenceDictionary --java-options "-Xmx${java_memory_param}G" \
-R ${output_file}

#####
#also need to adjust the variants file
#####

original_vcf=${emmer_wheat_modern_variants_original}
split_header=${extra_data_dir}/vcf_header_split.txt
split_vcf=${TMPDIR}/all_emmer_filtered_variants_header_to_SAMN04448013.split.vcf
split_vcf_edited_header=${TMPDIR}/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.vcf
split_vcf_edited_header_snp_only=${TMPDIR}/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf
output_vcf=${output_dir}/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

#make a temporary split vcf
cp $original_vcf ${split_vcf}.tmp.vcf

#for each chromosome, we take the ${split_vcf}.tmp.vcf file and adjust the chromosome names using awk and output to ${split_vcf}. Then we move the current $split_vcf file to ${split_vcf}.tmp.vcf and repeat the process for the next chromosome
while read line ; do
  chr=$(echo $line | awk '{print $1}')
  length=$(echo $line | awk '{print $2}')
  echo "editing $chr of length $length at $(date)"
  awk -F $'\t' '$1 == "'$chr'" && $2 <= 400000000 {$1 = "'$chr':1-400000000"} {OFS = FS} {print}' ${split_vcf}.tmp.vcf | \
  awk -F $'\t' '$1 == "'$chr'" && $2 >= 400000001 && $2 <= '$length' {$1 = "'$chr':400000001-'$length'";$2-=400000000} {OFS = FS} {print}' > $split_vcf
  mv $split_vcf ${split_vcf}.tmp.vcf
done < ${TMPDIR}/chr_names_lengths.txt

#finally, the ${split_vcf}.tmp.vcf can be moved to ${split_vcf}
mv ${split_vcf}.tmp.vcf ${split_vcf}

##### next we want to adjust the header and we will remove indels
echo "adding new header at $(date)"
#copy the new header to the new file location
cp $split_header $split_vcf_edited_header
#take all the contents of the split vcf after 37 lines (the old header size) and add onto the new header
tail -n+37 $split_vcf >> $split_vcf_edited_header

#use vcftools to remove indels
echo "removing indels at $(date)"
vcftools --vcf $split_vcf_edited_header \
--remove-indels \
--recode-INFO-all \
--recode \
--out $split_vcf_edited_header_snp_only

#copy back from TMPDIR
echo "copying split snp only vcf with edited header at $(date)"
cp $split_vcf_edited_header_snp_only* $output_vcf

#index feature file for GATK
gatk IndexFeatureFile --java-options "-Xmx${java_memory_param}G" \
--feature-file $output_vcf

#also zip and index this file for use with vcftools
bcftools convert ${output_vcf} --output ${output_vcf}.gz --output-type z
bcftools index ${output_vcf}.gz

#finally we will check coverage over these snps using a bedfile 
vcf2bed < ${output_vcf} | awk -v "OFS=\t" '{print $1, $2, $3}' > ${output_vcf}.bed
