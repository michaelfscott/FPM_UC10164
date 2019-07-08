args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==3) {
  stop("supply arguments: (1) extra_data_dir (2) traw_file (3) threshold_param")
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
#plot_export_prefix ="plots/"
threshold=as.numeric(args[3])/100

traw<-read_tsv(traw_file) %>% join_chromosome_halves_rename_samples() 
location<-load_metadata(metadata_filename)

###########
# fraction windows plot
###########

get_high_conc<-function(sample, threshold){
	MAF_by_window<-read_tsv(paste0("RData/",sample,"_window_50.tsv")) %>%
		tidyr::gather(key="Accession", value="MAF", -c(chromosome, physical_position))
	ungroup(MAF_by_window) %>% group_by(Accession) %>%
		summarise(Accession2=sample, high_conc=sum(MAF<threshold, na.rm=TRUE) / sum(!is.na(MAF)))
}
get_low_conc<-function(sample, threshold){
	MAF_by_window<-read_tsv(paste0("RData/",sample,"_window_50.tsv")) %>%
		tidyr::gather(key="Accession", value="MAF", -c(chromosome, physical_position))
	result<-ungroup(MAF_by_window) %>% group_by(chromosome, physical_position) %>%
		summarise(high_conc=sum(MAF<threshold, na.rm=TRUE)) %>%
		ungroup() %>% summarise(Accession=sample, unique_windows=sum(high_conc<=1), num_windows=n(), unique_windows_frac= unique_windows/num_windows) 
	return(result)
}

accession_cols<-location$Accession[location$Accession %in% colnames(traw) & ! location$Accession %in% exclude_samples]

all_high_conc<-lapply(accession_cols, get_high_conc, threshold=threshold)
high_conc_tbl<-purrr::reduce(all_high_conc, rbind)
all_low_conc<-lapply(accession_cols,get_low_conc, threshold=threshold)
low_conc_tbl<-purrr::reduce(all_low_conc, rbind)
write_tsv(low_conc_tbl, path=paste0("RData/low_conc_", threshold, ".tsv"))
write_tsv(high_conc_tbl, path=paste0("RData/high_conc_", threshold, ".tsv"))

}
