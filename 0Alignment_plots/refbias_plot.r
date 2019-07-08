library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

full_traw_path<-"../../5Variants/3Modern_Samples_Merge_bedfilter/mapQ20/merged.traw"
traw_full<-read_tsv(full_traw_path)
traw_full_geno<-traw_full[,7:ncol(traw_full)]
alts<-apply(traw_full_geno,2,sum,na.rm=TRUE)/(2*apply(!is.na(traw_full_geno),2,sum))
extension="all_sites"
tmp<-data.frame(name=names(alts), non_reference_allele_frac=alts) %>% tbl_df() %>%
	mutate(
		Accession=ifelse(grepl("Mt_", name)| grepl("Iraq_", name), paste0(word(name,1,sep="_"),"_", word(name,2,sep="_")), word(name,1,sep="_")), 
		type=ifelse(Accession %in% c("UC10164S1", "Outgroup"), Accession, "other"),
		extension=paste0(extension))

get_alts<-function(traw_path){
	traw<-read_tsv(traw_path)
	traw_geno<-traw[,7:ncol(traw)] %>% mutate(UC10164S2_UC10164S2=as.integer(UC10164S2_UC10164S2))
	apply(traw_geno,2,sum,na.rm=TRUE)/(2*apply(!is.na(traw_geno),2,sum))
}

get_data_frame<-function(mapQfilter, minDP_S1, extension){
	alts<-get_alts(paste0("../../5Variants/3Modern_Samples_Merge_",extension, "/mapQ", mapQfilter, "/S1minDP", minDP_S1, "/merged.traw"))
	plotter<-data.frame(name=names(alts), non_reference_allele_frac=alts) %>% tbl_df() %>%
		mutate(
			Accession=ifelse(grepl("Mt_", name)| grepl("Iraq_", name), paste0(word(name,1,sep="_"),"_", word(name,2,sep="_")), word(name,1,sep="_")), 
			type=ifelse(Accession %in% c("UC10164S1", "Outgroup"), Accession, "other"),
			extension=paste0(extension))
}

mapQfilters=c(20,25,30)
minDP_S1_levels=c(2,4)
for(i in seq(1:3)){
mapQfilter=mapQfilters[i]
for(j in seq(1:2)){
minDP_S1=minDP_S1_levels[j]

tmp1<-get_data_frame(mapQfilter, minDP_S1, "bedfilter")
tmp2<-get_data_frame(mapQfilter, minDP_S1, "bedfilter_no_indels")
tmp3<-get_data_frame(mapQfilter, minDP_S1, "bedfilter_no_indels_altbam")
tmp4<-get_data_frame(mapQfilter, minDP_S1, "bedfilter_altref")
tmp5<-get_data_frame(mapQfilter, minDP_S1, "bedfilter_no_indels_altref_altbam")

ext_from=c("all_sites","bedfilter", "bedfilter_no_indels", "bedfilter_no_indels_altbam", "bedfilter_altref", "bedfilter_no_indels_altref_altbam")
ext_to=c("1) all sites", paste0("2) aDNA depth >= ", minDP_S1), paste0("3) indel alignments removed, aDNA depth >= ", minDP_S1), paste0("5) alternate reads, indel alignments removed, aDNA depth >= ", minDP_S1), paste0("4) alternate reference genome, aDNA depth >= ", minDP_S1), paste0("6) alternate reads, alternate reference genome, indel alignments removed, aDNA depth >= ", minDP_S1))

refbias_data<-bind_rows(tmp,tmp1,tmp2,tmp3,tmp4,tmp5) %>%
	mutate(filter=plyr::mapvalues(extension, from=ext_from, to=ext_to))

#combine with metadata about location etc. 
location<-read.table("../../2Reference_Genomes/Zavitan_v2/metadata/origins.txt", stringsAsFactors=FALSE, header=TRUE) %>% tbl_df() %>% 
	rbind(list("UC10164S1","Egypt","UC10164S1", "dicoccum", "ancient_S1", "Ancient_Egyptian", 26.93, 31.48)) %>%
	rbind(list("UC10164S2","Egypt","UC10164S2", "dicoccum", "ancient_S2", "Ancient_Egyptian", 26.93, 31.48)) %>%
	rbind(list("Outgroup", NA, "Outgroup", "Outgroup", "Outgroup", "Outgroup", NA,NA))

plotter<-left_join(refbias_data, location)

customPalette <- c("#000000", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#000066", "#888888", "#56B4E9", "#999999")
label_names<- c("Ancient Egyptian", "Caucasus (Dom)", "Eastern Europe (Dom)", "Indian Ocean (Dom)", "Mediterranean (Dom)", "Northern Levant (Wild)", "Outgroup", "Southern Levant (Wild)")

refbias_plot_mapQ_region<-ggplot(plotter, aes(x= non_reference_allele_frac, fill=Region)) +
	geom_histogram(binwidth=0.005) +
	facet_wrap(~filter, ncol=1) +
	xlab("% non-reference alleles") +
	theme_minimal() +
	scale_fill_manual(labels= label_names , values= customPalette)+
	scale_color_manual(labels= label_names , values= customPalette) +
NULL

refbias_plot_mapQ_IS<-ggplot(plotter, aes(x= non_reference_allele_frac, fill=Improvement_status)) +
	geom_histogram(binwidth=0.005) +
	facet_wrap(~filter, ncol=1) +
	xlab("% non-reference alleles") +
	theme_minimal() +
	scale_fill_manual(values= customPalette[c(1,2,3,9,8)])+
NULL

ggsave(file=paste0("plots/refbias_plot_minDP", minDP_S1,"_mapQ",mapQfilter,".pdf"), refbias_plot_mapQ_IS, height=9, width=7)

#check there is actually a difference in altref:	
#filter(plotter, extension %in% c("bedfilter", "bedfilter_altref"), Accession=="UC10164S1")

}
}
