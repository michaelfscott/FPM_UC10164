#$ -cwd
#$ -N STEP1_convert_and_qpDstat_Southern_Levant
#$ -o STEP1_convert_and_qpDstat_Southern_Levant.log
#$ -l mem=1G
#$ -l h_rt=06:00:00 
#$ -S /bin/bash
#$ -j y
#$ -l tmpfs=5G
#$ -t 1-15

set -a

source ../1Input_scripts/input_data_locations.sh
source ../1Input_scripts/software.sh

reference=${project_dir}/2Reference_Genomes/Zavitan_v2_split/151210_zavitan_v2_pseudomolecules_split.fasta
known_variants=${project_dir}/2Reference_Genomes/Zavitan_v2_split/all_emmer_filtered_variants_header_to_SAMN04448013.split.edit.header.snp.vcf

#choose mapQ filter and samplename parameters across SGE_TASK_ID's
get_filters() {
  local TASK_ID=$1
  local mapQs="20 25 30"
  local filters=".bedfilter .bedfilter.noindels .bedfilter.posfilter.altref .bedfilter.noindels.posfilter.altBam .bedfilter.noindels.posfilter.altRef.altBam"
  local filter_extensions="bedfilter bedfilter_no_indels bedfilter_altref bedfilter_no_indels_altbam bedfilter_no_indels_altref_altbam"
  local mapQ_num=$((($TASK_ID+2)%3 + 1))
  mapQfilter=$(echo $mapQs | cut -d" " -f$mapQ_num)
  local filter_num=$((($TASK_ID+2)/3))
  filter=$(echo $filters | cut -d" " -f$filter_num)
  extension=$(echo $filter_extensions | cut -d" " -f$filter_num)
}

get_filters $SGE_TASK_ID

#for plink, we speficy the allele to count so that we count alternate alleles:
awk -v OFS='\t' '{print $1":"$2,$5}' ${known_variants} | \
grep -v "^#" > ${TMPDIR}/reference_alleles.txt
#we will remove the S2 sample from the genotypes used for AdmixTools analysis
remove_samples_file=${TMPDIR}/remove_samples_file.txt
echo "UC10164S2 UC10164S2" > $remove_samples_file

#individual and poplist files have been manually edited to give the relevant populations for testing
indivfile=indfiles/Split_Southern_Levant_indfile.ind
poplistfile=poplists/Split_Southern_Levant_poplist.txt

for minDP_S1 in 2 4; do

plink_input=${project_dir}/5Variants/3Modern_Samples_Merge_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/merged
plink_converted=${TMPDIR}/merged
convert_out_dir=${project_dir}/8AdmixTools/1Converted_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}
convert_out_prefix=${convert_out_dir}/merged
mkdir -p $convert_out_dir

plink \
--bfile $plink_input \
--remove $remove_samples_file \
--reference-allele  ${TMPDIR}/reference_alleles.txt \
--chr-set 30 \
--recode \
--out $plink_converted

convertf_parfile=${TMPDIR}/merged.par
echo "genotypename:    ${plink_converted}.ped" > $convertf_parfile
echo "snpname:         ${plink_converted}.map" >> $convertf_parfile
echo "indivname:       ${plink_converted}.ped" >> $convertf_parfile
echo "outputformat:    EIGENSTRAT" >> $convertf_parfile
echo "genotypeoutname: ${convert_out_prefix}.geno" >> $convertf_parfile
echo "snpoutname:      ${convert_out_prefix}.snp" >> $convertf_parfile
echo "indivoutname:    ${convert_out_prefix}.ind" >> $convertf_parfile
echo "familynames:     NO" >> $convertf_parfile
echo "noxdata:     YES" >> $convertf_parfile
echo "outputgroup:     YES" >> $convertf_parfile
echo "numchrom:     28" >> $convertf_parfile

convertf -p $convertf_parfile

#now perform qpDstat
prefix=$(basename $convert_out_prefix)
admixtools_out_dir=${project_dir}/8AdmixTools/2AdmixTools_out_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}
mkdir -p $admixtools_out_dir

#repeat across a variety of block sizes
blgsizes="0.5 0.1 0.05 0.01 0.005 0.001"
blgsize_extensions="50mb 10mb 5mb 1mb 500kb 100kb"

for i in {1..6}; do 

blgsize=$(echo $blgsizes | cut -d" " -f$i )
blgsize_extension=$( echo $blgsize_extensions | cut -d" " -f$i)

qpDstat_parfile=${TMPDIR}/qpDstat.${blgsize_extension}.par
echo "DIR:   ${convert_out_dir}" > $qpDstat_parfile
echo "SSS:   $prefix" >> $qpDstat_parfile
echo "indivname:    $indivfile" >> $qpDstat_parfile
echo "snpname:      DIR/SSS.snp" >> $qpDstat_parfile
echo "genotypename: DIR/SSS.geno" >> $qpDstat_parfile
echo "popfilename:  $poplistfile" >> $qpDstat_parfile
echo "printsd:  YES" >> $qpDstat_parfile
echo "blgsize:  $blgsize" >> $qpDstat_parfile
echo "numchrom:   28" >> $qpDstat_parfile

qpDstat -p $qpDstat_parfile > ${admixtools_out_dir}/Split_Southern_Levant_qpDstat.${blgsize_extension}.results.txt

done

done
