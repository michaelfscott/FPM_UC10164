args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==3) {
  stop("supply arguments: (1) extra_data_dir (2) traw_file (3) export_file")
} else {

source("../0Rfunctions/functions.r")

###########
# load data
###########

#extra_data_dir="../1Input_scripts/extra_data"
extra_data_dir=args[1]
metadata_filename=paste0(extra_data_dir, "/origins_with_location.txt")
blast_res<-read.table(paste0(extra_data_dir, "/array_blast_result7.txt"))
#traw_file="../../5Variants/3Modern_Samples_Merge_bedfilter_no_indels_altref_altbam/mapQ30/S1minDP2/merged.traw"
traw_file=args[2]
#export_file=randomization_test.txt
export_file=args[3]

traw<-read_tsv(traw_file) %>% join_chromosome_halves_rename_samples()
location<-load_metadata(metadata_filename)
outliers<-load_diversity_scan_data(paste0(extra_data_dir,"/"))

###########
# concordance within outlier regions
###########

sample="UC10164S1"
#window size in diversity scan
window_size<-2000000
increments_per_window=2

window_positions<-get_window_positions(traw, window_size, increments_per_window)

#round all the positions and remove columns for joining
traw_tbl<-tbl_df(traw) %>% mutate(rounded=plyr::round_any(position, window_size/2, f=ceiling)) %>% select(-c("CHR", "SNP", "COUNTED", "ALT"))
#join together datasets
outliers_genos<-full_join(select(outliers, -c("CHR", "Chr")), 
		full_join(window_positions , traw_tbl, by=c("chromosome", "rounded")), 
	by=c("chromosome", "Start_pos")) %>% 
	select(-"Start_pos") %>%
	arrange(chromosome, position) 

compare<-compare_with_target(outliers_genos, sample)

concordances<-concordance_by_outlier_type(compare)

###########
# diversity scan outliers significance
###########

num_windows<-nrow(window_positions)
num_replicates=1000

#perform randomization test
set.seed(123456789)
result<-lapply(1:num_replicates, function(x)(random_test_statistic(compare, num_windows)))
export<-do.call(bind_rows, result)
write_tsv(export, export_file)
#export<-read_tsv("randomization_test.txt")

exclude_samples=c("UC10164S2", "Outgroup")
empirical<-concordances %>%
	filter(!Outlier_type=="other_loci", !Accession %in% exclude_samples) %>%
	left_join(location) %>%
	group_by(Outlier_type) %>% 
	summarise(wild_dom_diff=
		mean(concordance[Improvement_status=="domesticated"], na.rm=TRUE) -
		mean(concordance[Improvement_status=="wild"], na.rm=TRUE))
write_tsv(empirical, paste0(export_file, ".empirical.txt"))

group_by(export, Outlier_type) %>%
	summarise(
		FST_significance=(sum(wild_dom_diff>pull(empirical[1,2]))+1)/(num_replicates+1),
		PiRaio_significance=(sum(wild_dom_diff>pull(empirical[2,2])+1)/(num_replicates+1)),
		TajD_significance=(sum(wild_dom_diff>pull(empirical[3,2]))+1)/(num_replicates+1))

}