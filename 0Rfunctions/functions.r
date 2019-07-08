library(ggplot2)
library(stringr)
library(reshape2)
library(readr)
library(grid)
library(gtable)
library(gridExtra)
library(dplyr)
library(ape)

#colour palette
customPalette <- c("#000000", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#000066", "#56B4E9")
Accession_plot_order <- c("PI467008","Tabigha","PI470962","PI466950","PI471060","PI471021","PI538673","Gamla","PI471038","PI428122","PI466946","Qazerin","PI466957","Nesher","PI428132","PI352322","Mt_Hermon","Mt_Gerizim","PI487264","Iraq_41","PI428129","PI428066","PI428084","PI538666","PI428054","PI428088","PI428057","PI538642","PI428072","PI428070","PI428036","PI428025","PI538631","PI254158","PI254169","PI254180","PI470739","PI470738","PI94661","PI470737","PI326312","PI94640","PI606325","PI352329","PI94741","PI377658","PI264964","PI434995","PI352347","PI355454","PI355496","PI352357","PI352348","PI352367","PI352352","PI182743","PI355498","PI352361","PI191091","PI276007","PI532302","PI322232","PI319868","PI319869","UC10164S1")
exclude_samples=c("UC10164S2", "Outgroup")

replace_na_with_median<-function(vector){
vector[is.na(vector)]<-median(vector, na.rm=TRUE)
return(vector)
}

get_chr_startpos<-function(traw){
  start_positions_per_chr<-traw %>% group_by(CHR) %>%
    summarize(
    startpos= as.numeric(word(word(SNP[1],2,sep=":"),1,sep="-"))-1,
    endpos= as.numeric(word(word(SNP[1],2,sep=":"),2,sep="-")),
    chromosome=word(SNP[1],1,sep=":"))
  return(start_positions_per_chr)
}

join_chromosome_halves_rename_samples<-function(traw){
  start_positions_per_chr<-get_chr_startpos(traw) 
  samplenames<-colnames(traw)[7:ncol(traw)]
  edited_samplenames<-ifelse(grepl("Mt_", samplenames)| grepl("Iraq_", samplenames), paste0(word(samplenames,1,sep="_"),"_", word(samplenames,2,sep="_")), word(samplenames,1,sep="_"))
  traw_new<-left_join(traw,start_positions_per_chr, by="CHR") %>% ungroup() %>%
    mutate(position=POS+startpos)
  #re-order columns
  column_order<-c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", samplenames)
  traw_new <-traw_new[,column_order]
  colnames(traw_new)<-c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", edited_samplenames)
  return(traw_new)
}

load_diversity_scan_data<-function(data_dir_diversity){
	Fst_outliers<-read.csv(paste0(data_dir_diversity, "FST_outliers.csv"), header=TRUE, stringsAsFactors=FALSE)
	Fst_outliers$Chr<-paste0("chr",word(Fst_outliers$Chr,1,sep="_"))
	PiRatio_outliers<-read.csv(paste0(data_dir_diversity, "PiRatio_outliers.csv"), header=TRUE, stringsAsFactors=FALSE)
	PiRatio_outliers$Chr<-paste0("chr",word(PiRatio_outliers$Chr,1,sep="_"))
	TajD_outliers<-read.csv(paste0(data_dir_diversity, "Dom_TajD_outliers.csv"), header=TRUE, stringsAsFactors=FALSE)
	TajD_outliers$Chr<-paste0("chr",word(TajD_outliers$Chr,1,sep="_"))
	overlap_outliers<-read.csv(paste0(data_dir_diversity, "Overlap_among_diversity_scans.csv"), header=TRUE, stringsAsFactors=FALSE)
	overlap_outliers$Chr<-paste0("chr",word(overlap_outliers$Chr,1,sep="_"))
	tbl_df(full_join(full_join(Fst_outliers, PiRatio_outliers, by=c("Chr", "Start_pos")), TajD_outliers, by=c("Chr", "Start_pos"))) %>% mutate(CHR=as.numeric(as.factor(Chr)), chromosome=Chr)
}

load_metadata<-function(metadata_filename){
	read_delim(metadata_filename, delim=" ") %>%
	rbind(list("UC10164S1","Egypt","UC10164S1", "dicoccum", "domesticated", "Ancient_Egyptian", 26.93, 31.48)) %>%
	rbind(list("UC10164S2","Egypt","UC10164S2", "dicoccum", "domesticated", "Ancient_Egyptian", 26.93, 31.48)) %>%
	rbind(list("Outgroup", NA, "Outgroup", "Outgroup", "Outgroup", "Outgroup", NA,NA)) %>%
	mutate(
		#relabel "Region" of two samples that seem to fall into the wrong group genetically and geographically
		Region=ifelse(Label %in% c("WE-19", "WE-20"), "Southern_Levant", Region),
		#create "type" of both improvment status and region
		type=paste(Improvement_status, Region, sep="_")
	)
}

get_window_positions<-function(traw, window_size, increments_per_window){
  start_positions_per_chr<-get_chr_startpos(traw)
  chromosome_lengths<-group_by(start_positions_per_chr, chromosome) %>% summarise(end_CHR=max(endpos)) 
  chromosome_lengths[chromosome_lengths$chromosome=="Un",2]<-480980714
  window_positions<-chromosome_lengths %>%
    rowwise() %>%
    do(
      data.frame(
        chromosome = .$chromosome, 
        Start_pos = seq(1000, .$end_CHR+ 1000, by = window_size/2),#(seq(0, .$end_CHR, by = window_size/increments_per_window) + 1000), 
        rounded = seq(window_size/increments_per_window, .$end_CHR+ window_size/increments_per_window, by = window_size/increments_per_window)
      )) %>% ungroup()
  return(window_positions)
}

compare_with_target<-function(outliers_genos, target_samplename){
	tidyr::gather(outliers_genos, key="Accession", value="genotype", -c("FST", "Dom_to_Wild_Pi_Ratio", "Dom_TajD", "position", target_samplename, "chromosome", "rounded")) %>%
	group_by(Accession, chromosome, rounded) %>%
	summarise(
		FST=!is.na(FST[1]),
		Dom_to_Wild_Pi_Ratio =!is.na(Dom_to_Wild_Pi_Ratio[1]),
		Dom_TajD =!is.na(Dom_TajD[1]),
		matching=sum(UC10164S1== genotype, na.rm=TRUE), 
		num_calls=sum(!is.na(UC10164S1) & !is.na(genotype)))
}

concordance_by_outlier_type_window<-function(compare){
	#in this part I use lead to check if the next window is also an outlier. If so, then only the focal window_size/2 region is counted. Else, both that window and the next window are counted. This is because we are going to sum up across all the windows so we don't want to count twice. 
	ungroup(compare) %>% 
	group_by(Accession, chromosome) %>% 
	arrange(Accession, chromosome, rounded) %>% #must be in the right order
	mutate(
		FST_matching=ifelse(!FST,NA,ifelse(lead(FST), matching, matching+lead(matching))),
		FST_num_calls=ifelse(!FST,NA,ifelse(lead(FST), num_calls, num_calls +lead(num_calls))),
		Dom_to_Wild_Pi_Ratio_matching=ifelse(! Dom_to_Wild_Pi_Ratio,NA,ifelse(lead(Dom_to_Wild_Pi_Ratio), matching, matching+lead(matching))),
		Dom_to_Wild_Pi_Ratio_num_calls=ifelse(! Dom_to_Wild_Pi_Ratio,NA,ifelse(lead(Dom_to_Wild_Pi_Ratio), num_calls, num_calls +lead(num_calls))),
		Dom_TajD_matching=ifelse(! Dom_TajD,NA,ifelse(lead(Dom_TajD), matching, matching+lead(matching))),
		Dom_TajD_num_calls=ifelse(! Dom_TajD,NA,ifelse(lead(Dom_TajD), num_calls, num_calls +lead(num_calls))),
		other_matching=ifelse(FST | Dom_to_Wild_Pi_Ratio | Dom_TajD,NA,ifelse(lag(FST | Dom_to_Wild_Pi_Ratio | Dom_TajD), NA, matching)),
		other_num_calls=ifelse(FST | Dom_to_Wild_Pi_Ratio | Dom_TajD,NA,ifelse(lag(FST | Dom_to_Wild_Pi_Ratio | Dom_TajD), NA, num_calls)))
}

concordance_by_outlier_type<-function(compare){
	conc<-concordance_by_outlier_type_window(compare)
	ungroup(conc) %>% group_by(Accession) %>%
summarise(FST=sum(FST_matching[FST], na.rm=TRUE)/sum(FST_num_calls[FST], na.rm=TRUE),
	PiRatio=sum(Dom_to_Wild_Pi_Ratio_matching[Dom_to_Wild_Pi_Ratio], na.rm=TRUE)/sum(Dom_to_Wild_Pi_Ratio_num_calls[Dom_to_Wild_Pi_Ratio], na.rm=TRUE),
	TajD=sum(Dom_TajD_matching[Dom_TajD], na.rm=TRUE)/sum(Dom_TajD_num_calls[Dom_TajD], na.rm=TRUE),
	other_loci=sum(other_matching, na.rm=TRUE)/sum(other_num_calls, na.rm=TRUE)) %>% 
	tidyr::gather(key="Outlier_type", value="concordance", -c("Accession"))
}

dom_minor_allele<-function(traw_regions, location){
	accession_cols<-location$Accession[location$Accession %in% colnames(traw_regions)]
	domesticated_columns<-accession_cols[(location$Improvement_status=="domesticated" & ! location$Region=="Ancient_Egyptian")]
	switch_minor_allele<-apply(traw_regions[,domesticated_columns],1,sum,na.rm=TRUE)/apply(!is.na(traw_regions[,domesticated_columns]),1,sum)>1
	#convert all scores 
	traw_dom_MAF<-traw_regions
	traw_dom_MAF[switch_minor_allele, accession_cols]<-abs(traw_regions[switch_minor_allele, accession_cols]-2)
	return(traw_dom_MAF)
}

MAF_by_SNP_window<-function(traw_dom_MAF, snp_window_size, window_move_size, location){
	accession_cols<-location$Accession[location$Accession %in% colnames(traw_dom_MAF)]
	tidyr::gather(traw_dom_MAF, key="Accession", value="genotype", accession_cols) %>%
		  group_by(chromosome, Accession) %>%
		arrange(chromosome, Accession, position) %>%
		mutate(
		num_calls_roll=zoo::rollapply(!is.na(genotype), snp_window_size, sum, by= window_move_size, fill=NA, na.rm=TRUE),
		minor_alleles_roll=zoo::rollapply(genotype, snp_window_size, sum, by= window_move_size, fill=NA, na.rm=TRUE),
		MAF= minor_alleles_roll/(2*num_calls_roll),
		physical_position = zoo::rollapply(position , snp_window_size, median, by= window_move_size, fill=NA, na.rm=TRUE)/1000000
    ) %>%
  filter(!is.na(physical_position))
} 

similarity_by_SNP_window<-function(traw, sample, snp_window_size, window_move_size, location){
	#convert all scores to count minor allele against target sample
	switch_minor_allele<-as.logical(traw[,sample]==2)
	heterozygous_genotype<-as.logical(traw[,sample]==1)
	missing_genotype<-is.na(switch_minor_allele)
	switch_minor_allele[missing_genotype]<-FALSE
	heterozygous_genotype[missing_genotype]<-FALSE
	traw_MAF<-traw
	traw_MAF[switch_minor_allele, accession_cols]<-abs(traw[switch_minor_allele, accession_cols]-2)
	traw_MAF[missing_genotype, accession_cols]<-NA
	traw_MAF[heterozygous_genotype, accession_cols]<-NA
	accession_cols<-location$Accession[location$Accession %in% colnames(traw_MAF)]
	tidyr::gather(traw_MAF, key="Accession", value="genotype", accession_cols) %>%
		  group_by(chromosome, Accession) %>%
		arrange(chromosome, Accession, position) %>%
		mutate(
		num_calls_roll=zoo::rollapply(!is.na(genotype), snp_window_size, sum, by= window_move_size, fill=NA, na.rm=TRUE),
		minor_alleles_roll=zoo::rollapply(genotype, snp_window_size, sum, by= window_move_size, fill=NA, na.rm=TRUE),
		MAF= minor_alleles_roll/(2*num_calls_roll),
		physical_position = zoo::rollapply(position , snp_window_size, median, by= window_move_size, fill=NA, na.rm=TRUE)/1000000
    ) %>%
  filter(!is.na(physical_position))
} 

MAF_by_physical_window<-function(traw_dom_MAF, window_size, increments_per_window, location){
	window_positions<-get_window_positions(traw_dom_MAF, window_size, increments_per_window)
	accession_cols<-location$Accession[location$Accession %in% colnames(traw_dom_MAF)]
	mutate(traw_dom_MAF , rounded=plyr::round_any(position, window_size/2, f=ceiling)) %>%
	full_join(window_positions , by=c("chromosome", "rounded")) %>%
	tidyr::gather(key="Accession", value="genotype", accession_cols) %>%
		group_by(chromosome, Accession, rounded) %>%
		summarise(num_calls=sum(!is.na(genotype)), 
			minor_alleles=sum(genotype, na.rm=TRUE)) %>%
		arrange(chromosome, Accession, rounded) %>%
		mutate(
			num_calls_roll=zoo::rollapply(num_calls, 2, sum, by= 1, fill=NA, na.rm=TRUE),
			minor_alleles_roll=zoo::rollapply(minor_alleles, 2, sum, by= 1, fill=NA, na.rm=TRUE),
			MAF= minor_alleles_roll/(2*num_calls_roll),
			physical_position = zoo::rollapply(rounded ,2, min, by=1, fill=NA, na.rm=TRUE)/1000000
		) %>%
  filter(!is.na(physical_position))
} 

max_MAF_by_region<-function(MAF_by_physical, columns, chromosome_name, lower_bound_mb, upper_bound_mb, threshold){
	filter(MAF_by_physical, Accession %in% columns) %>%
		group_by(chromosome, rounded) %>%
		summarize(max_dom_MAF=max(MAF, na.rm=TRUE)) %>%
		filter(max_dom_MAF<threshold) %>%
		filter(chromosome== chromosome_name, rounded/1000000 > lower_bound_mb,  rounded/1000000 <upper_bound_mb) 
}

shift_compare<-function(compare, variable, shift){
	length<-nrow(compare)
	compare_tmp<-compare
	compare_tmp[,variable]<-lag(pull(compare,variable), shift)
	compare_tmp[1:shift, variable]<-compare[(length-shift+1):length, variable]
	return(compare_tmp)
}

shift_and_get_test_statistic<-function(compare, shift_FST, shift_PiRatio, shift_TajD, exclude_samples=c("UC10164S2", "Outgroup")){
	#append the shift statistic at the end
	shifts<-data.frame(Outlier_type=c("FST", "PiRatio", "TajD"), shift=c(shift_FST, shift_PiRatio, shift_TajD))
	#now perform shift
	shift_compare(compare, "FST", shift_FST) %>% 
	shift_compare( "Dom_to_Wild_Pi_Ratio", shift_PiRatio) %>%
	shift_compare( "Dom_TajD", shift_TajD) %>%
	concordance_by_outlier_type() %>%
	filter(!Outlier_type=="other_loci", !Accession %in% exclude_samples) %>%
	left_join(location) %>%
	group_by(Outlier_type) %>% 
	summarise(wild_dom_diff=
		mean(concordance[Improvement_status=="domesticated"], na.rm=TRUE) -
		mean(concordance[Improvement_status=="wild"], na.rm=TRUE)) %>%
	full_join(shifts)
}

random_test_statistic<-function(compare, num_windows, exclude_samples=c("UC10164S2", "Outgroup")){
	shift_FST=sample(1:num_windows,1)
	shift_PiRatio =sample(1:num_windows,1)
	shift_TajD=sample(1:num_windows,1)
	shift_and_get_test_statistic(compare, shift_FST, shift_PiRatio, shift_TajD, c("UC10164S2", "Outgroup"))
}

grab_grob <- function(){
  grid.echo()
  grid.grab()
}
