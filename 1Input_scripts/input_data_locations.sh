#all data will be stored in the project dir
project_dir=/lustre/projects/MAGIC_WHEAT/UC10164/
#the fastqs for UC10164 should be found within the following directory
raw_fastq_dir=${project_dir}/1raw_fastqs/
#the location of the psuedochromosome fasta of the Emmer wheat reference genome should be specified here
emmer_wheat_reference_original=${project_dir}/2Reference_Genomes/Zavitan_v2/151210_zavitan_v2_pseudomolecules.fasta
#the location of the vcf of emmer wheat variants should be specified here
emmer_wheat_modern_variants_original=${project_dir}/2Reference_Genomes/Zavitan_v2/all_emmer_filtered_variants_header_to_SAMN04448013.vcf
#the location of the psuedochromosome fasta of the bread wheat reference genome should be specified here
bread_wheat_reference_original=${project_dir}/2Reference_Genomes/IWGSC_refseq_v1/161010_Chinese_Spring_v1.0_pseudomolecules.fasta
#a vcf for the outgroup 
outgroup_vcf=${project_dir}/1outgroup_vcf/Outgroup_AB_genome.vcf

#extra input data (e.g., metadata for modern samples) is included in a subdirectory specified here
extra_data_dir=${project_dir}/scripts/1Input_scripts/extra_data/
