library(readr)
library(dplyr)
library(stringr)
library(reshape2)
library(ggplot2)

##### coverage statistics

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
	chromosome_subgenome<-read_tsv(paste0(S1_path, "bedcov/chromosomal.mapQ", 0,".bedcov"), col_names=columnNames) %>%
	mutate(subgenome=str_replace(str_extract(chromosome, "chr[0-9][A,B,D]"), "chr[0-9]","")) %>% select(chromosome, subgenome)
	S1<-list()
	S2<-list()
	for(i in seq(1,7)){
		min_mapQ<-min_mapQs[i]
		S1[[i]]<-get_coverage_chromosomal(path=paste0(S1_path, "bedcov/chromosomal.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
		S2[[i]]<-get_coverage_chromosomal(path=paste0(S2_path, "bedcov/chromosomal.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
	}
	names(S1)<-min_mapQs
	names(S2)<-min_mapQs
	bind_rows(list(S1=bind_rows(S1, .id = "min_mapQ") , S2=bind_rows(S2, .id = "min_mapQ")), .id="sample")
}

get_all_coverage_SNP<-function(S1_path, S2_path, columnNames, min_mapQs, by_subgenome=FALSE){
	chromosome_subgenome<-read_tsv(paste0(S1_path, "bedcov/chromosomal.mapQ", 0,".bedcov"), col_names=columnNames) %>%
	mutate(subgenome=str_replace(str_extract(chromosome, "chr[0-9][A,B,D]"), "chr[0-9]","")) %>% select(chromosome, subgenome)
	S1<-list()
	S2<-list()
	for(i in seq(1,7)){
		min_mapQ<-min_mapQs[i]
		S1[[i]]<-get_coverage_SNP(path=paste0(S1_path, "bedcov/SNP.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
		S2[[i]]<-get_coverage_SNP(path=paste0(S2_path, "bedcov/SNP.mapQ", min_mapQ,".bedcov"), columnNames=columnNames, chromosome_subgenome, by_subgenome)
	}
	names(S1)<-min_mapQs
	names(S2)<-min_mapQs
	bind_rows(list(S1=bind_rows(S1, .id = "min_mapQ") , S2=bind_rows(S2, .id = "min_mapQ")), .id="sample")
}


##### misincorporation_patterns_

get_misincorporation<-function(S1_path, S2_path){
	tmp<-list()
	paths=c(S1_path, S2_path)
	for(i in seq(1,2)){
		tmp[[i]]<-read.table(paste0(paths[i],"map_damage/misincorporation.txt"), header=TRUE)
	}
	names(tmp)<-c("S1", "S2")
	full<-bind_rows(tmp, .id="sample")
	transitions<-tbl_df(full) %>% 
	  group_by(End,Pos, sample) %>% 
	  summarise(GtoA=sum(G.A)/sum(G), CtoT=sum(C.T)/sum(C)) %>%
	  data.frame()
	mis<-melt(transitions, id.vars=c("End", "Pos", "sample"))
	colnames(mis)<-c("End","Position","sample", "Misincorporation","Frequency")
	mis[mis$End=="5p","Position"]<-(-mis[mis$End=="5p","Position"])
	mis$End<-gsub("3p","3 prime end", gsub("5p","5 prime end", mis$End))
	mis$Misincorporation<-gsub("GtoA","G to A", gsub("CtoT","C to T", mis$Misincorporation))
	return(mis)
}

##### fragment sizes

get_fragment_lengths<-function(S1_path, S2_path){
	tmp<-list()
	paths=c(S1_path, S2_path)
	for(i in seq(1,2)){
		tmp[[i]]<-read.table(paste0(paths[i],"stats/readlengths.txt"), col.names=c("Length", "Occurences"))
	}
	names(tmp)<-c("S1", "S2")
	lengths<-bind_rows(tmp, .id="sample") %>% 
	  group_by(sample, Length) %>% 
	  summarise(Occurences=sum(Occurences)) %>%
	  data.frame()
	return(lengths)
}


######### PLOT FUNCTIONS #########

plot_lengths<-function(lengths){
	ggplot(lengths, aes(Length, Occurences/1000000)) + 
	  geom_col(width=1) +
	  facet_grid(sample~., scales = "free_y") +
	  theme(legend.position="top") +
	  xlab("Fragment Length") +
	  ylab("Frequency (millions)") +
	  ggtitle("Fragment Sizes") +
	  theme_minimal() + 
	  theme(
	    panel.grid.minor = element_blank(),
	    plot.title = element_text(hjust = 0.5),
	    strip.background = element_blank(),
		strip.placement = "outside",
		panel.border = element_rect(colour = "black", fill=NA)) +
	  NULL
}

plot_mis<-function(mis){
	ggplot(mis, aes(x= Position, y=Frequency, colour=Misincorporation)) + 
	  geom_line() + 
	  facet_grid(sample~End, scales = "free_x") +ylim(0,0.1) +
	  xlab("distance from fragment end") +
	  ylab("frequency") +
	  ggtitle("Misincorporations versus Emmer Reference") +
	  theme_minimal() + 
	  theme(
	    legend.position="top",
	    legend.title = element_blank(),
	    plot.title = element_text(hjust = 0.5),
	    panel.grid.minor = element_blank(),
	    strip.background = element_blank(),
		strip.placement = "outside",
		panel.border = element_rect(colour = "black", fill=NA)) +
	  NULL  
}

plot_bread_subgenome<-function(tbl_bread){
	ggplot(tbl_bread, aes(x= min_mapQ, y=coverage, fill=min_mapQ)) +
	  geom_bar(stat = "identity", position = position_dodge()) +
	  scale_fill_grey(start=0.8, end=0) +
	  facet_grid(sample~subgenome, scale="free_y") +
	  ggtitle("Alignment to Hexaploid Bread Reference") +
	  theme_minimal() + 
	  xlab("minimum mapping quality") +
	  theme(
	    legend.position = "none",
	    plot.title = element_text(hjust = 0.5),
	    panel.grid.minor = element_blank(),
	    panel.grid.major.x = element_blank(),
	    strip.background = element_blank(),
		strip.placement = "outside",
		panel.border = element_rect(colour = "black", fill=NA)) +
	  NULL
}

plot_emmer_full_SNP<-function(tbl_emmer){
	mutate(tbl_emmer, location=ifelse(location=="SNP_sites", "exonic SNP sites", "all sites")) %>% #change for labelling
	ggplot(aes(x= min_mapQ, y=coverage, fill=min_mapQ)) +
	  geom_bar(stat = "identity", position = position_dodge()) +
	  scale_fill_grey(start=0.8, end=0) +
	  facet_grid(sample~ location, scale="free_y") +
	  ggtitle("Alignment to Emmer Reference") +
	  theme_minimal() + 
	  xlab("minimum mapping quality") +
	  theme(
	    plot.title = element_text(hjust = 0.5),
	    legend.position = "none",
	    panel.grid.minor = element_blank(),
	    panel.grid.major.x = element_blank(),
	    strip.background = element_blank(),
		strip.placement = "outside",
		panel.border = element_rect(colour = "black", fill=NA)) +
	  NULL
}