<h2> Ancient Egyptian Emmer Wheat </h2>

These scripts are for replication of the analyses in Scott et al. (2019) "Whole genome sequence from 3,000-year-old Egyptian emmer wheat reveals dispersal and domestication history". 

<h3>Overview </h3>

<p>Scripts are organised into groups that should be run sequentially. For example, 2Reference_Genomes is used to prepare the reference genomes, which is necessary for Alignment and Variant Calling (3Alignments_And_Variant_Calls).
</p>
<p>Within each directory, each script begins with STEP[N]. All scripts with the same number [N] are independent of one another. However, STEP[N] may require that all STEP[N-1] jobs have been completed. That is, these steps should be run sequentially.
</p>

<h3>COMPUTING ENVIRONMENT </h3>

<p>These scripts were tested on the UCL Research Computing HPC cluster, which uses a Sun Grid Engine scheduler. The flags for submission via SGE are at the head of each submission script.
</p>

<p>Using SGE, parallel jobs can be submitted using the flag #$ -t 1-n, which submits n jobs, each of which has an environmental variable "$SGE_TASK_ID". Therefore, to replicate these scripts in a different environment, they may need to be edited accordingly (e.g., adding a for loop of SGE_TASK_ID from {1..n}. 
</p>

<p>
Another environmental variable is provided when submitting scripts: $TMPDIR. This specifies a temporary local directory in which intermediate files can be stored. On this HPC system, $TMPDIR is automatically deleted upon job completion. In other environments, TMPDIR may need to be manually specified and/or removed. 
</p>

<p>The required software is specified/loaded in "1Input_scripts/software.sh". This will have to be modified to reflect the local install locations.
</p>

<h3>INPUT DATA </h3>

<p>Sequence data for two samples are deposited at the ENA under the sample objects ERS3135152 and ERSXXXXXXXXX
</p>

<p>The Wild Emmer Wheat (T. turgidum ssp.dicoccoides) reference genome psuedochromosomes (Avni et al. 2017), 151210_zavitan_v2_pseudomolecules.fasta, can be downloaded from https://wewseq.wixsite.com/consortium. 

<p>In addition, variants from exome capture of a panel of wild and domesticated emmer wheats (Avni et al. 2017), all_emmer_filtered_variants_header_to_SAMN04448013.vcf, are available from the same consortium web address.
</p>

<p>The SAM format specifications impose an upper limit on the length attribute allowable for chromosomes. It is therefore common practice to split wheat chromosomes. In our analysis, we split each of the psuedochromosomes at 400mb. The psuedochromosomes can be split using 2Reference_Genomes_scripts/STEP1_Split_Emmer_Reference.sh.
</p>


<p>To check for contamination with modern bread wheat, we align against the IWGSC Triticum aestivum refseq1 reference genome. pseudochromosomes can be obtained from ensembl (https://plants.ensembl.org/Triticum_aestivum/Info/Index) or URGI (https://urgi.versailles.inra.fr/download/iwgsc/IWGSC_RefSeq_Assemblies/v1.0/)
</p>

<p>The locations into which these files have been downloaded can be specified in the "1Input_scripts/input_data_locations.sh". 
</p>


<h3> VARIANTS </h3>

<p>To skip the alignment and variant calling steps. The variants can be downloaded from
</p>

<p>This vcf also includes genotypes called from diploid outgroup species using data accession numbers SRR4010671 and SRR4010672 (Triticum urartu) and SAMEA2342530 (Aegilops speltoides) used for A and B subgenome outgroup genotypes, respectively. 
</p>
