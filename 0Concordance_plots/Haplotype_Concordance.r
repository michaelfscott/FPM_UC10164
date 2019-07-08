args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==4) {
  stop("supply arguments: (1) extra_data_dir (2) traw_file (3) export_dir (4) sample_number")
} else {

source("../0Rfunctions/functions.r")

###########
# load data
###########

#extra_data_dir="../1Input_scripts/extra_data"
extra_data_dir=args[1]
metadata_filename=paste0(extra_data_dir, "/origins_with_location.txt")
#traw_file="../../5Variants/3Modern_Samples_Merge_bedfilter_no_indels_altref_altbam/mapQ30/S1minDP2/merged.traw"
traw_file=args[2]
#export_dir="RData"
export_dir=args[3]
#i=1
i=as.numeric(args[4])

traw<-read_tsv(traw_file) %>% join_chromosome_halves_rename_samples() 
location<-load_metadata(metadata_filename)

###########
# similarity plots for sample
###########

snp_window_size=50
window_move_size=1
exclude_samples=c("UC10164S2", "Outgroup")
accession_cols<-location$Accession[location$Accession %in% colnames(traw) & !location$Accession %in% exclude_samples]
traw<-select(traw, -exclude_samples)

#for(i in 1:length(accession_cols)) {
	sample<-accession_cols[i]
	
	MAF_by_window<-similarity_by_SNP_window(traw, sample, snp_window_size, window_move_size, location)
	
	output<-tidyr::spread(select(ungroup(MAF_by_window), c(chromosome, physical_position, Accession, MAF)), key=Accession, value=MAF)
	output_path=paste0(export_dir,"/",sample,"_window_", snp_window_size, ".tsv")
	write_tsv(output, output_path)
#}

}
