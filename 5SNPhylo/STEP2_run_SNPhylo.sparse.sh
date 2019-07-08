#parameters
num_bootstraps=1000
LD_threshold=0.1

project_dir=/home/ucbtmf1/MOUNT_MAGIC_WHEAT/UC10164/

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

get_filters $1

for minDP_S1 in 2 4; do

input_vcf=${project_dir}/6SNPhylo/1converted_vcf_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/SNPhylo_incl_UC10164S1.vcf
input_vcf2=${project_dir}/6SNPhylo/1converted_vcf_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/SNPhylo_excl_UC10164S1.vcf
outdir=${project_dir}/6SNPhylo/2SNPhylo_result_${extension}/mapQ${mapQfilter}/S1minDP${minDP_S1}/
mkdir -p $outdir

cd $outdir

${snphylo} -v $input_vcf \
-c $minDP_S1 \
-m 0.1 \
-M 0.1 \
-l $LD_threshold \
-a 28 \
-o Outgroup \
-b \
-B $num_bootstraps 
#-m is the MAF threshold
#-M is the missing rate

${snphylo} -v $input_vcf2 \
-c $minDP_S1 \
-m 0.1 \
-M 0.1 \
-l $LD_threshold \
-a 28 \
-o Outgroup \
-b \
-B $num_bootstraps 
#-m is the MAF threshold
#-M is the missing rate

done

#######
#now exclude UC10164S1 entirely
#######

input_vcf=${project_dir}/6SNPhylo/1converted_vcf_${extension}/mapQ${mapQfilter}/excludeS1/SNPhylo_excl_UC10164S1.vcf
outdir=${project_dir}/6SNPhylo/2SNPhylo_result_${extension}/mapQ${mapQfilter}/excludeS1/
mkdir -p $outdir

cd $outdir

${snphylo} -v $input_vcf \
-c $minDP_S1 \
-m 0.1 \
-M 0.1 \
-l $LD_threshold \
-a 28 \
-o Outgroup \
-b \
-B $num_bootstraps
#-m is the MAF threshold
#-M is the missing rate

