library(readr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)

get_coverage_chromosomal<-function(path, columnNames, chromosome_subgenome, by_subgenome=FALSE){
	if(by_subgenome){
		tbl<-read_tsv(path, col_names=columnNames) %>%
		full_join(chromosome_subgenome, by="chromosome") %>%
		mutate(subgenome=ifelse(is.na(subgenome), "Un", subgenome)) %>%
		group_by(subgenome) %>% 
		summarise(raw_mapped=sum(raw/1000000), dedup_mapped=sum(dedup/1000000), length=sum((end/1000000)), raw_coverage= raw_mapped/length, dedup_coverage= dedup_mapped/length) 
		out<-data.frame(
		status=c(rep("raw", nrow(tbl)), rep("dedup", nrow(tbl))), 
		subgenome=rep(tbl$subgenome, 2),
		coverage=c(tbl$raw_coverage, tbl$dedup_coverage))
		return(out)
	} else {
		tbl<-read_tsv(path, col_names=columnNames) %>%
		summarise(raw_mapped=sum(raw/1000000), dedup_mapped=sum(dedup/1000000), length=sum((end/1000000)), raw_coverage= raw_mapped/length, dedup_coverage= dedup_mapped/length) 
		out<-data.frame(
		status=c(rep("raw", nrow(tbl)), rep("dedup", nrow(tbl))), 
		coverage=c(tbl$raw_coverage, tbl$dedup_coverage))
		return(out)
	}
}

get_coverage_SNP<-function(path, columnNames, chromosome_subgenome, by_subgenome=FALSE){
	if(by_subgenome){
		tbl<-read_tsv(path, col_names=columnNames) %>%
		left_join(chromosome_subgenome, by="chromosome") %>%
		mutate(subgenome=ifelse(is.na(subgenome), "Un", subgenome)) %>%
		group_by(subgenome) %>% 
		summarise(raw_coverage= sum(as.numeric(raw))/n(), dedup_coverage= sum(as.numeric(dedup))/n())
		out<-data.frame(
		status=c(rep("raw", nrow(tbl)), rep("dedup", nrow(tbl))), 
		subgenome=rep(tbl$subgenome, 2),
		coverage=c(tbl$raw_coverage, tbl$dedup_coverage))
		return(out)
	} else {
		tbl<-read_tsv(path, col_names=columnNames) %>%
		summarise(raw_coverage= sum(as.numeric(raw))/n(), dedup_coverage= sum(as.numeric(dedup))/n())
		out<-data.frame(
		status=c(rep("raw", nrow(tbl)), rep("dedup", nrow(tbl))), 
		coverage=c(tbl$raw_coverage, tbl$dedup_coverage))
		return(out)
	}
}

get_all_coverage_chromosomal<-function(S1_path, S2_path, columnNames, min_mapQs, by_subgenome=FALSE){
	chromosome_subgenome<-read_tsv(paste0(S1_path, "chromosomal.mapQ", 0,".bedcov"), col_names=columnNames) %>%
	mutate(subgenome=str_replace(str_extract(chromosome, "chr[0-9][A,B,D]"), "chr[0-9]","")) %>% select(chromosome, subgenome)
	S1<-list()
	S2<-list()
	for(i in seq(1,7)){
		min_mapQ<-min_mapQs[i]
		S1[[i]]<-get_coverage_chromosomal(path=paste0(S1_path, "chromosomal.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
		S2[[i]]<-get_coverage_chromosomal(path=paste0(S2_path, "chromosomal.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
	}
	names(S1)<-min_mapQs
	names(S2)<-min_mapQs
	bind_rows(list(S1=bind_rows(S1, .id = "min_mapQ") , S2=bind_rows(S2, .id = "min_mapQ")), .id="sample")
}

get_all_coverage_SNP<-function(S1_path, S2_path, columnNames, min_mapQs, by_subgenome=FALSE){
	chromosome_subgenome<-read_tsv(paste0(S1_path, "chromosomal.mapQ", 0,".bedcov"), col_names=columnNames) %>%
	mutate(subgenome=str_replace(str_extract(chromosome, "chr[0-9][A,B,D]"), "chr[0-9]","")) %>% select(chromosome, subgenome)
	S1<-list()
	S2<-list()
	for(i in seq(1,7)){
		min_mapQ<-min_mapQs[i]
		S1[[i]]<-get_coverage_SNP(path=paste0(S1_path, "SNP.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
		S2[[i]]<-get_coverage_SNP(path=paste0(S2_path, "SNP.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
	}
	names(S1)<-min_mapQs
	names(S2)<-min_mapQs
	bind_rows(list(S1=bind_rows(S1, .id = "min_mapQ") , S2=bind_rows(S2, .id = "min_mapQ")), .id="sample")
}

S1_path="../4Alignments/1Align_Bread_Collapsed/Merged/UC10164S1/bedcov/"
S2_path="../4Alignments/1Align_Bread_Collapsed/Merged/UC10164S2/bedcov/"
min_mapQs<-c(0,1,10,20,25,30,35)
columnNames<-c("chromosome", "start", "end", "raw", "dedup")
tbl<-get_all_coverage_chromosomal(S1_path, S2_path, columnNames, min_mapQs, by_subgenome=TRUE)

filter(tbl, status=="raw", !subgenome=="Un") %>%
mutate(status_and_mapQ=paste0(status , min_mapQ)) %>%
ggplot(aes(x= min_mapQ, y=coverage, alpha=min_mapQ)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_alpha_discrete(range=c(0.3,1)) +
  facet_grid(sample~subgenome, scale="free_y") +
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
	strip.placement = "outside",
	panel.border = element_rect(colour = "black", fill=NA)) +
  NULL

S1_path="../4Alignments/1Align_Emmer_Collapsed_Trimmed/Merged/UC10164S1/bedcov/"
S2_path="../4Alignments/1Align_Emmer_Collapsed_Trimmed/Merged/UC10164S2/bedcov/"
min_mapQs<-c(0,1,10,20,25,30,35)
tbl_emmer<-bind_rows(list(
	SNP_sites= get_all_coverage_SNP(S1_path, S2_path, columnNames, min_mapQs) , 
	all= get_all_coverage_chromosomal(S1_path, S2_path, columnNames, min_mapQs)), .id="location")

filter(tbl_emmer, status=="raw") %>%
ggplot(aes(x= min_mapQ, y=coverage, alpha=min_mapQ)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  scale_alpha_discrete(range=c(0.3,1)) +
  facet_grid(sample~ location, scale="free_y") +
  theme_minimal() + 
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    strip.background = element_blank(),
	strip.placement = "outside",
	panel.border = element_rect(colour = "black", fill=NA)) +
  NULL
