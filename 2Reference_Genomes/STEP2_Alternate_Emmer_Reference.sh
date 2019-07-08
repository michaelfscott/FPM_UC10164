#job name
#$ -N STEP2_Alternate_Emmer_Reference
#$ -o STEP2_Alternate_Emmer_Reference.log
# Set the working directory 
#$ -cwd
#memory and runtime options 
#$ -l mem=40G
#$ -l h_rt=48:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=90G

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

java_memory_param=10

reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
alternate_reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.alternate.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

#make bed file from split_vcf
#awk -v FS='\t' -v OFS='\t' '/^[^#]/ { $2 = $2 - 1; print $1, $2, $2 }' $known_variants > ${TMPDIR}/split_SNPs.bed
#try without using zero-based bed file
awk -v FS='\t' -v OFS='\t' '/^[^#]/ { print $1, $2, $2 }' $known_variants > ${TMPDIR}/split_SNPs.bed
#also requires a tped file
#first convert vcf to plink using vcftools
#for plink, we speficy the allele to count so that we count alternate alleles:
awk -v OFS='\t' '{print $1":"$2,$5}' ${known_variants} | \
grep -v "^#" > ${TMPDIR}/reference_alleles.txt
#should also specify the chromosome map so that the chromosomes are numbered sequentially:
awk -v OFS='\t' '{print $1, NR}' ${reference}.fai > ${TMPDIR}/chromosome_map.txt
#output the tped file
vcftools --vcf $known_variants \
--plink --chrom-map ${TMPDIR}/chromosome_map.txt \
--out ${TMPDIR}/split_vcf_plink
plink \
--file ${TMPDIR}/split_vcf_plink \
--recode transpose \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--out ${TMPDIR}/split_vcf_plink

#python3 requred for refbias scripts (2.7 loaded for vcftools and plink).
module load python2/recommended
python2 ${REFBIAS_PATH}/modify_refseq_different_allele.py ${TMPDIR}/split_vcf_plink.tped $reference ${TMPDIR}/split_SNPs.bed > $alternate_reference

#index new reference genome
echo "bwa indexing reference at $(date)"
bwa index ${alternate_reference}
echo "samtools faidx reference at $(date)"
samtools faidx ${alternate_reference}
echo "GATK indexing reference at $(date)"
gatk CreateSequenceDictionary --java-options "-Xmx${java_memory_param}G" -R ${alternate_reference}

