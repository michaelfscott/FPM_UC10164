#here we specify various standard bioinformatic tools for analysis of NGS data
module load bwa/0.7.12/gnu-4.9.2
module load samtools/1.3.1/gnu-4.9.2
module load java/1.8.0_92 gatk/4.0.8.0
module load bedtools/2.25.0
module load perl/5.22.0
module load vcftools/0.1.15/gnu-4.9.2
module load plink/1.90b3.40

#refbias script location:
REFBIAS_PATH=/lustre/projects/MAGIC_WHEAT/UC10164/scripts/1Input_scripts/refbias/
#SNPhylo script location:
SNPhylo_PATH=/home/ucbtmf1/packages/snphylo/SNPhylo/
#ADMIXTURE
PATH=/home/ucbtmf1/packages/admixture_linux-1.3.0/:$PATH
#AdmixTools
PATH=/home/ucbtmf1/packages/AdmixTools/bin/:$PATH
#add htslib (tabix and bgzip)
PATH=/home/ucbtmf1/packages/htslib/:$PATH
#add bcftools to path
PATH=/home/ucbtmf1/packages/bcftools/:$PATH
#we also use the fastp toolkit to split paired-end fastq files
PATH=/home/ucbtmf1/packages/fastp_v0.19.7/:$PATH
#we also add the adapter removal program to the path
PATH=/home/ucbtmf1/packages/adapterremoval-2.2.2/build/:$PATH
#we also use the fastx toolkit
PATH=/home/ucbtmf1/packages/fastx_toolkit_0.0.13_binaries_Linux_2.6_amd64/:$PATH
#the mapdamage function is in ~/.local/bin. The dependencies (e.g., pysam, R) must also be installed
PATH=/home/ucbtmf1/.python2local/bin/:$PATH
#vcf2bed is included in
PATH=/home/ucbtmf1/packages/bin/:$PATH


#To load r on this system, we must also modify the compilers
module unload compilers
module unload mpi
module load r/recommended

