#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if (!length(args)==4) {
  stop("supply arguments: (1) metadata_filename (2) input_fam_file (3) admixture_prefix (4) export_dir")
} else {
	
source("../0Rfunctions/functions.r")

#metadata_filename="../../2Reference_Genomes/Zavitan_v2/metadata/origins.txt"
metadata_filename=args[1]
#input_fam_file="../../7Admixture/1LD_pruned_bedfilter_no_indels_altref_altbam/mapQ30/S1minDP2/merged.pruned.fam"
input_fam_file=args[2]
#admixture_prefix="../../7Admixture/2ADMIXTURE_out_bedfilter_no_indels_altref_altbam/mapQ30/S1minDP2/merged"
admixture_prefix=args[3]
#export_dir="plots"
export_dir=args[4]

order<-read_delim(input_fam_file, col_names=FALSE, delim=" ") %>% mutate(Accession=X1) %>% select(Accession) # read in order of .fam file to match ADMIXTURE output with location data. 
location<-load_metadata(metadata_filename) %>% right_join(order)

#classify "Tabigha_group" so that the output plot can be properly coloured
location_tmp<-location
location_tmp$Region[location_tmp$Accession %in% c("Tabigha", "PI467008", "PI470962", "PI466950", "PI471021")]<-"Tabigha_group"

#replace column names so that they are consistent across K and reflect most common 'region' 
replace_col_name<-function(tbl, region, target_name, location_tmp){
	tmp<-names(which.max(apply(tbl[location_tmp$Region== region,],2,sum)))
	tmp2<-colnames(tbl)== target_name
	colnames(tbl)[colnames(tbl)==tmp]<-target_name
	colnames(tbl)[tmp2]<-tmp
	return(tbl)
}
convert_all_col_names<-function(tbl,location_tmp){
	tbl<-replace_col_name(tbl,"Indian_Ocean", "V1", location_tmp)
	tbl<-replace_col_name(tbl,"Tabigha_group", "V7", location_tmp)
	tbl<-replace_col_name(tbl,"Southern_Levant", "V6", location_tmp)
	tbl<-replace_col_name(tbl,"Northern_Levant", "V5", location_tmp)
	tbl<-replace_col_name(tbl,"Mediterranean", "V4", location_tmp)
	tbl<-replace_col_name(tbl,"Caucasus", "V3", location_tmp)
	tbl<-replace_col_name(tbl,"Eastern_Europe", "V2", location_tmp)
}
#read in data
tbl2<-read.table(paste0(admixture_prefix,"/merged.pruned.2.Q"))
tbl3<-read.table(paste0(admixture_prefix, "/merged.pruned.3.Q"))
tbl4<-read.table(paste0(admixture_prefix, "/merged.pruned.4.Q"))
tbl5<-read.table(paste0(admixture_prefix, "/merged.pruned.5.Q"))
tbl6<-read.table(paste0(admixture_prefix, "/merged.pruned.6.Q"))
tbl7<-read.table(paste0(admixture_prefix, "/merged.pruned.7.Q"))
#convert column names
tbl2<-convert_all_col_names(tbl2,location_tmp)
tbl3<-convert_all_col_names(tbl3,location_tmp)
tbl4<-convert_all_col_names(tbl4,location_tmp)
tbl5<-convert_all_col_names(tbl5,location_tmp)
tbl6<-convert_all_col_names(tbl6,location_tmp)
tbl7<-convert_all_col_names(tbl7,location_tmp)

#melt the data
df<-data.frame(Accession= location$Accession, Improvement_status= location$Improvement_status, Region= location$Region, K2=tbl2,K3=tbl3, K4=tbl4, K5=tbl5, K6=tbl6, K7=tbl7)
df_melted<-melt(df, id.vars = c("Accession","Improvement_status", "Region"))
# first part is the K value used, second part is the population
df_melted$Kvalue=word(df_melted$variable,1,sep="\\.")
df_melted$variable=word(df_melted$variable,2,sep="\\.")

#for colouring the labels
colours<-sapply(as.numeric(as.factor(df$Region[order(df$Region)])), function(x)(customPalette[x]))

#to get appropriate labels
location$label_number<-as.numeric(gsub("UC10164S1","10164",gsub("[A-Z][A-Z]\\-","",location$Label)))
plot_order<-order(factor(location$Region, levels=rev(c("Ancient_Egyptian", "Indian_Ocean", "Mediterranean", "Eastern_Europe", "Caucasus", "Northern_Levant", "Southern_Levant"))), location$label_number)
plot_order<-plot_order[c(11:14,10,1:9,15:65)]
#location$Accession[plot_order]
colours<-sapply(as.numeric(as.factor(location $Region)), function(x)(customPalette[x]))[plot_order]

df_melted$Kvalue<-gsub("K", "K=",df_melted$Kvalue)

output<-ggplot(df_melted, aes(x= factor(Accession, levels=df$Accession[plot_order]), y=value, fill=variable)) + 
  geom_bar(stat = "identity") +
  ylab("Ancestry") +
  xlab("Accession") +
  theme_minimal() +
  theme(legend.position='none',
    axis.ticks.x= element_line(),
    axis.text.x = element_blank(),
    axis.title.x= element_text(),
    axis.title.y= element_text(),
    axis.text.y = element_text(colour = colours)) +
  coord_flip() +
  scale_fill_manual(values = c(customPalette[-1][c(3,2,1,4,5,6)],"#CCFFFF")) +
  facet_wrap(~ Kvalue, nrow=1) +
  NULL

ggsave(paste0(export_dir,"/Admixture.pdf"), output, height=7, width=7)

cross_v<-read_tsv(paste0(admixture_prefix,"/CV_error.txt"), col_names=c("CV_error"))
cross_v$K_param=seq(2,7) 
cross_v<-select(cross_v, c("K_param", "CV_error"))
write_csv(data.frame(cross_v), paste0(export_dir,"/CV_error.csv"))
#ggplot(cross_v, aes(x=K_param,y=cross_v[,2])) + geom_point()

}
