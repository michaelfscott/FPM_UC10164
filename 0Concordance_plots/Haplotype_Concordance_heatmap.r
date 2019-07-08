args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==3) {
  stop("supply arguments: (1) extra_data_dir (2) traw_file (3) plot_export_prefix")
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
plot_export_prefix =args[3]

traw<-read_tsv(traw_file) %>% join_chromosome_halves_rename_samples() 
location<-load_metadata(metadata_filename)

###########
# fraction windows plot
###########

sample="UC10164S1"
sample="Tabigha"

MAF_by_window<-read_tsv(paste0("RData/",sample,"_window_50.tsv")) %>%
	tidyr::gather(key="Accession", value="MAF", -c(chromosome, physical_position))

threshold=0.05
ungroup(MAF_by_window) %>% group_by(Accession) %>%
	summarise(sum(MAF<threshold, na.rm=TRUE) / sum(!is.na(MAF))) %>% print(n=100)

min_MAF<-ungroup(MAF_by_window) %>% group_by(chromosome, position) %>%
	summarise(min_MAF=min(MAF,na.rm=TRUE)) 

filter(MAF_by_window, Accession=="UC10164S1")

all_high_conc<-lapply(location$Accession[1:2],get_high_conc, threshold=threshold)
high_conc_matrix<-purrr::reduce(all_high_conc, right_join, by = "Accession")

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

for(threshold in c(seq(0.01,0.1,by=0.01), seq(0.2,0.5,by=0.1))){
	all_high_conc<-lapply(accession_cols, get_high_conc, threshold=threshold)
	high_conc_tbl<-purrr::reduce(all_high_conc, rbind)
	all_low_conc<-lapply(accession_cols,get_low_conc, threshold=threshold)
	low_conc_tbl<-purrr::reduce(all_low_conc, rbind)
	write_tsv(low_conc_tbl, path=paste0("RData/low_conc_", threshold, ".tsv"))
	write_tsv(high_conc_tbl, path=paste0("RData/high_conc_", threshold, ".tsv"))
}s

#all_high_conc<-lapply(location$Accession,get_high_conc)
all_high_conc<-lapply(location$Accession[1:2], get_high_conc, threshold=threshold)
high_conc_tbl<-purrr::reduce(all_high_conc, rbind)
get_low_conc<-function(sample, threshold){
	MAF_by_window<-read_tsv(paste0("RData/",sample,"_window_50.tsv")) %>%
		tidyr::gather(key="Accession", value="MAF", -c(chromosome, physical_position))
	result<-ungroup(MAF_by_window) %>% group_by(chromosome, physical_position) %>%
		summarise(high_conc=sum(MAF<threshold, na.rm=TRUE)) %>%
		ungroup() %>% summarise(Accession=sample, unique_windows=sum(high_conc<=1), num_windows=n(), unique_windows_frac= unique_windows/num_windows) 
	return(result)
}
all_low_conc<-lapply(location$Accession,get_low_conc, threshold=threshold)
low_conc_tbl<-purrr::reduce(all_low_conc, rbind)

#save(all_low_conc, all_high_conc, file=paste0("RData/conc_",threshold,".RData"))
#}
#for(threshold in c(seq(0.01,0.1,by=0.01), seq(0.2,0.5,by=0.1))){
#load(file= paste0("RData/conc_",threshold,".RData"))
#print(which(names(sort(all_low_conc, decreasing=TRUE))=="UC10164"))
#print(which(names(sort(all_low_conc[location$Improvement_status=="domesticated"], decreasing=TRUE))=="UC10164"))
#}

threshold=0.05
load(file= paste0("RData/conc_",threshold,".RData"))
load(file= paste0("RData/ncol_res.RData"))

high_conc_matrix<-reduce(all_high_conc, right_join, by = "Accession")
high_conc_matrix<-data.frame(high_conc_matrix)
rownames(high_conc_matrix)<-high_conc_matrix$Accession
high_conc_matrix<-high_conc_matrix[,-1] #remove accession column
colnames(high_conc_matrix)<-rownames(high_conc_matrix)[1:2]
concMelted<-melt(as.matrix(high_conc_matrix))

location$label_number<-as.numeric(gsub("UC10164","0",gsub("[A-Z][A-Z]\\-","",location$Label)))
plot_order<-order(factor(location$Region, levels=rev(c("Ancient_Egyptian", "Indian_Ocean", "Mediterranean", "Eastern_Europe", "Caucasus", "Northern_Levant", "Southern_Levant"))), location$label_number)
plot_order<-plot_order[c(11:14,10,1:9,15:65)]
colours<-sapply(as.numeric(as.factor(location $Region)), function(x)(customPalette[x]))[plot_order]

p1<-ggplot(concMelted,aes(
  x=factor(Var1, levels= levels(concMelted$Var1)[rev(plot_order)]),
  y=factor(Var2, levels= levels(concMelted$Var2)[plot_order]),
  fill=value)) +
  geom_tile() +
  labs(x="",y="Accession")+
  scale_x_discrete(expand=c(0,0))+
  scale_fill_gradientn(name="Fraction 'concordant' sliding windows", limits=c(0,1),
						colours = c("blue","red","yellow","white"), 
						values = c(0,0.33,0.66,1),
						guide= guide_colourbar(title.position="top", title.hjust=0.5)) +
  theme(legend.position="bottom", 
        legend.key.width=unit(.1,"npc"),legend.key.height=unit(.02,"npc"), 
        axis.text.y = element_text(colour = colours),
        axis.text.x = element_blank())


tmp<-data.frame(location, all_low_conc= all_low_conc)
p2<-ggplot(tmp, aes(x= all_low_conc/ncol_res), 
  y=factor(tmp$Region, levels=rev(c("Ancient_Egyptian", "Indian_Ocean", "Mediterranean", "Eastern_Europe", "Caucasus", "Northern_Levant", "Southern_Levant"))), 
  colour=Region,
  pch=Improvement_status)) +
  geom_jitter(height=0.1) +
  scale_colour_manual(values = customPalette) +
  scale_shape_manual(values=c(16, 15))+
  theme(legend.position="none", 
    axis.title.y = element_blank()) + 
  labs(x=paste0("Fraction sliding windows that are \'unique\'")) +
#  labs(x=paste0("Fraction of 50-snp windows with no comparisons >",100*(1-threshold),"% concordant"),y="Region") +
  NULL  

pdf(file=paste0(plots_chr_scan_dir,"high_concordance_windows_TH",threshold,".pdf"), height=8, width=8)
p1
dev.off()
  
pdf(file=paste0(plots_chr_scan_dir,"unique_windows_TH",threshold,".pdf"), height=2.5, width=5)
p2
dev.off()

library(plyr)
concMelted2<- concMelted
concMelted2$Var2<-mapvalues(concMelted2$Var2, from=location$Accession, to=location$Region)

p3<-concMelted2 %>% filter(value<1) %>%
  filter(Var2 %in% c("Northern_Levant", "Southern_Levant")) %>%
  filter(Var1 %in% location[location$Improvement_status=="domesticated","Accession"]) %>%
ggplot(aes(
  y=value, 
  x=factor(Var1, levels= rev(levels(concMelted$Var1)[plot_order])), 
  colour=factor(Var2, levels=c("Ancient_Egyptian", "Indian_Ocean", "Mediterranean", "Eastern_Europe", "Caucasus", "Northern_Levant", "Southern_Levant"))) ) +
  geom_jitter(width=0.1) +
  stat_summary(fun.y = mean, geom = "errorbar", 
             aes(ymax = ..y.., ymin = ..y.., group = Var2),
             width = 0.8, linetype = "solid") +
  scale_colour_manual(values = customPalette[c(6,7)]) +
  scale_shape_manual(values=c(16, 15)) +
  theme(legend.position="none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(colour = rev(colours)[1:sum(location$Improvement_status=="domesticated")], angle = 90, hjust=1)) + 
  labs(y=paste0("Fraction 'similar' sliding windows")) +
  NULL  

pdf(file=paste0(plots_chr_scan_dir,"high_concordance_windows_dom_vs_wild_TH",threshold,".pdf"), height=3, width=9)
p3
dev.off()



################################
######## Concordance with UC10164 split by chromosome
################################

#UPDATE ONCE PLOT_ORDER IS DECIDED

#total_concordance<-(ungroup(conc) %>% group_by(variable) %>% summarise(concordance=sum(matching)/sum(num_calls)) %>% data.frame())[,2]

#export1<-conc %>% summarise(concordance=sum(matching)/sum(num_calls)) %>% ungroup() %>%
#  tidyr::spread(key=CHR, value=concordance) 

#location2<-location[export1[,"Accession"],]
#location2$label_number<-as.numeric(gsub("UC10164","10164",gsub("[A-Z][A-Z]\\-","",location2$Label)))
#plot_order<-order(factor(location2$Region, levels=rev(c("Ancient_Egyptian", "Indian_Ocean", "Mediterranean", "Eastern_Europe", "Caucasus", "Northern_Levant", "Southern_Levant"))), location2$label_number)
#plot_order<-plot_order[c(11:14,10,1:9,15:65)]

#export2<-data.frame(export1[rev(plot_order[-65]),])

#write.csv(export2, file=paste0(plot_export_prefix,"concordance_by_chr.csv"),row.names=FALSE, quote=FALSE)

