<h2> Ancient Egyptian Emmer Wheat </h2>

These scripts are for replication of the analyses in Scott et al. (2019) "Whole genome sequence from 3,000-year-old Egyptian emmer wheat reveals dispersal and domestication history". 


<h3>COMPUTING ENVIRONMENT </h3>
#
#These scripts were run on the UCL Research Computing HPC cluster, which uses a Sun Grid Engine scheduler. The flags for submission via SGE are at the head of each submission script.
#
#Using SGE, parallel jobs can be submitted using the flag #$ -t 1-n, which submits n jobs, each of which has an environmental variable "$SGE_TASK_ID. Therefore, to replicate these scripts in a different environment, they may need to be edited accordingly (e.g., adding a for loop of SGE_TASK_ID from {1..n}. 
#
#Another environmental variable is provided when submitting scripts: $TMPDIR. This specifies a temporary local directory in which intermediate files can be stored. On this HPC system, $TMPDIR is automatically deleted upon job completion. In other environments, TMPDIR may need to be manually specified and/or removed. 
#
#The required software is specified/loaded in "1Input_scripts/software.sh". This will have to be modified to reflect the local install locations. 
#
##### INPUT DATA #####
#
#Sequence data for two samples are deposited at the ENA under the sample objects ERS3135152 and ERSXXXXXXXXX
#
#The Wild Emmer Wheat (T. turgidum ssp.dicoccoides) reference genome psuedochromosomes (Avni et al. 2017), 151210_zavitan_v2_pseudomolecules.fasta, can be downloaded from https://wewseq.wixsite.com/consortium. 
#In addition, variants from exome capture of a panel of wild and domesticated emmer wheats (Avni et al. 2017), all_emmer_filtered_variants_header_to_SAMN04448013.vcf, are available from the same consortium web address.
#The SAM format specifications impose an upper limit on the length attribute allowable for chromosomes. It is therefore common practice to split wheat chromosomes. In our analysis, we split each of the psuedochromosomes at 400mb. 
#
#We also use data from SRR4010671 and SRR4010672 (Triticum urartu) and SAMEA2342530 (Aegilops speltoides) to get an outgroup genotype for the A and B subgenomes, respectively. We provide the .vcf files obtained from alignment of this sequence data to the (split) emmer reference genome.
#
#To check for contamination with modern bread wheat, we align against the IWGSC Triticum aestivum refseq1 reference genome. pseudochromosomes can be obtained from ensembl (https://plants.ensembl.org/Triticum_aestivum/Info/Index) or URGI (https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/)
#
#
#The locations into which these files have been downloaded can be specified in the "1Input_scripts/input_data_locations.sh" script
#
#
##### VARIANTS #####
#
#To skip the alignment and variant calling steps. The variants can be downloaded from
#
#
