args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==4) {
  stop("supply arguments: (1) extra_data_dir (2) traw_file (3) traw_regions_file (4) plot_export_prefix")
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
#traw_regions_file="../../5Variants/3Modern_Samples_Merge_bedfilter_no_indels_altref_altbam/mapQ30/merged.traw"
traw_regions_file=args[3]
#plot_export_prefix ="plots/"
plot_export_prefix =args[4]

traw<-read_tsv(traw_file) %>% join_chromosome_halves_rename_samples()
traw_regions<-read_tsv(traw_regions_file) %>% join_chromosome_halves_rename_samples()
location<-load_metadata(metadata_filename)
outliers<-load_diversity_scan_data(paste0(extra_data_dir,"/"))

###########
#plot parameters
###########

#outlier_plot_parameters
outlier_plot_colours<-customPalette[c(3,7)]
outlier_plot_pch<-c(16, 18)
outlier_plot_alpha=0.65
outlier_plot_theme<- theme_minimal() +
  	theme(legend.position="none", 
  		axis.title.x=element_blank(),
  		panel.grid.major.x =element_blank(),
  		strip.placement = "outside") 
outlier_ylab= "concordance with UC10164"
#to allow facets to be named using expressions
facet_names <- list(
  'FST'=expression(F[ST]*" outliers"),
  'PiRatio'=expression(pi[D]/pi[W]*" outliers"),
  'TajD'="Tajima's D outliers",
  'other_loci'="other loci"
)
facet_labeller <- function(variable,value){
  return(facet_names[value])
}

#scan plot parameters
padding_lower=30
padding_upper=30
scan_ylim=c(-0.05,1)
scan_size=0.3
scan_plot_theme<- theme_minimal() +
  	theme(legend.position="top", 
  		legend.title=element_blank()) 
scan_xlab<-"median physical position of window (mb)"
scan_ylab<-"frequency of minor allele"
scan_plot_colours<-customPalette[c(1,3,7)]
scan_plot_alpha<-c(1, 1, 0.25)

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

#filter_out accessions that are genotyped but not required in this comparison
outlier_plotter<-filter(concordances, !Accession %in% c("UC10164S2", "Outgroup")) %>%
	left_join(location)
#set order of facets for plotting
outlier_plotter$Outlier_type<-factor(outlier_plotter$Outlier_type, levels=c("FST", "PiRatio", "TajD", "other_loci"))

#produce plot
outlier_plots<-tbl_df(outlier_plotter) %>%
ggplot(aes(y= concordance, x=Improvement_status, colour= Improvement_status, shape=Improvement_status)) + 
  geom_jitter(alpha= outlier_plot_alpha, width=0.2) +
  facet_wrap(~ Outlier_type, ncol=4, labeller= facet_labeller) + 
  scale_colour_manual(values = outlier_plot_colours) +
  scale_shape_manual(values= outlier_plot_pch) +
  ylab(outlier_ylab) +
  outlier_plot_theme +
  NULL

###########
# similarity plots for UC10164
###########

########### positions of interest

Br_loci<-data.frame(marker_name=c("TtBtr1-A\n","TtBtr1-B\n"),physical_position=c(61639327/1000000, 97628063/1000000), chromosome=c("chr3A","chr3B"))
chr3A_Br<-data.frame(marker_name=c("TtBtr1-A\n"),physical_position=c(61639327/1000000), chromosome=c("chr3A"))
chr3B_Br<-data.frame(marker_name=c("TtBtr1-B\n"),physical_position=c(97628063/1000000), chromosome=c("chr3B"))
chr4B_physical_position<-apply(blast_res[blast_res$V2=="chr4B",9:10],1,mean)
chr4B_marker_name<-gsub("IWB","",blast_res[blast_res$V2=="chr4B",1])
chr4B_blast<-data.frame(marker_name=chr4B_marker_name,physical_position= chr4B_physical_position/1000000, chromosome=rep("chr4B",length(chr4B_physical_position)))
#remove those outside the plot range
chr4B_blast<-chr4B_blast[chr4B_blast$physical_position<580 & chr4B_blast$physical_position>420,]
chr4B_QTL_rect<-data.frame(xmin=min(chr4B_blast$physical_position), xmax=max(chr4B_blast$physical_position), ymin=0, ymax=1)
chr4B_QTL_region<-data.frame(x=c(475.3916, 540.1921), y=c(1,1))

################################
################ scan
################################

#traw_regions includes all the calls for UC10164S1. traw only includes those with at least depth 2. 
#the following mutates the calls for UC10164S1 in traw_regions to include only those with at least depth 2
traw_sample<-traw[,c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT", "UC10164S1")] %>% 
	mutate(UC10164S1_DP2=UC10164S1) %>% select(-UC10164S1)
traw_regions_DP2<-full_join(traw_regions, traw_sample, by= c("chromosome", "CHR", "SNP", "position", "COUNTED", "ALT")) %>%
	mutate(UC10164S1=UC10164S1_DP2) %>% select(-UC10164S1_DP2)

#make sure the minor allele in the domesticated population is labelled 0
accession_cols<-colnames(traw_regions_DP2)[colnames(traw_regions_DP2)  %in% location$Accession]
domesticated_columns<-colnames(traw_regions_DP2)[colnames(traw_regions_DP2)  %in% location$Accession[(location$Improvement_status=="domesticated" & ! location$Region=="Ancient_Egyptian")]]
switch_minor_allele<-apply(traw_regions_DP2[,domesticated_columns],1,sum,na.rm=TRUE)/apply(!is.na(traw_regions_DP2[,domesticated_columns]),1,sum)>1
#convert all scores 
traw_dom_MAF<-traw_regions_DP2
traw_dom_MAF[switch_minor_allele, accession_cols]<-abs(traw_regions_DP2[switch_minor_allele, accession_cols]-2)

#using SNP windows
snp_window_size=150
window_move_size=10

MAF_by_SNP<-MAF_by_SNP_window(traw_dom_MAF, snp_window_size, window_move_size, location)

to_plot<-tbl_df(location) %>% right_join(MAF_by_SNP) %>%
	mutate(Improvement_status=ifelse(Accession==sample, "ancient", Improvement_status)) %>%
	filter(!Accession %in% c("UC10164S2", "Outgroup"))

chr3A_plot<-tbl_df(to_plot) %>% 
	filter(chromosome=="chr3A", physical_position > chr3A_Br$physical_position-padding_lower,  physical_position <chr3A_Br$physical_position+ padding_upper) %>%
ggplot(aes(x = physical_position, y = MAF)) + 
	geom_vline(data= chr3A_Br, aes(xintercept= physical_position), size=0.75) +
	geom_text(data= chr3A_Br, aes(x=physical_position, label=marker_name, y=0.85), angle=90, size=3) +
	geom_line(aes(color = Improvement_status, group = Accession, alpha= Improvement_status), size=scan_size) +
	facet_grid(chromosome~., scales="free") +
	ylim(scan_ylim) +
	xlab(scan_xlab) +
	ylab(scan_ylab) +
	scale_color_manual(values = scan_plot_colours) +
	scale_alpha_manual(values = scan_plot_alpha) +
	scan_plot_theme +
	NULL

chr3B_plot<-tbl_df(to_plot) %>% 
	filter(chromosome=="chr3B", physical_position > chr3B_Br$physical_position-padding_lower,  physical_position <chr3B_Br$physical_position+ padding_upper) %>%
ggplot(aes(x = physical_position, y = MAF)) + 
	geom_vline(data= chr3B_Br, aes(xintercept= physical_position), size=0.75) +
	geom_text(data= chr3B_Br, aes(x=physical_position, label=marker_name, y=0.85), angle=90, size=3) +
	geom_line(aes(color = Improvement_status, group = Accession, alpha= Improvement_status), size=scan_size) +
	facet_grid(chromosome~., scales="free") +
	ylim(scan_ylim) +
	xlab(scan_xlab) +
	ylab(scan_ylab) +
	scale_color_manual(values = scan_plot_colours) +
	scale_alpha_manual(values = scan_plot_alpha) +
	scan_plot_theme +
	NULL
  
chr4B_plot_tbl<-filter(to_plot, chromosome == "chr4B", physical_position > 515-50,  physical_position <515+ 50)  

chr4B_plot<-ggplot() + 
	geom_line(data= chr4B_QTL_region, aes(x=x, y=y), color="black", size=0.75) +
	geom_point(data= chr4B_blast, aes(x=physical_position, y=1), color="black", shape="|", size=6) +
	geom_text(data= chr4B_blast, aes(x=physical_position, label=marker_name, y=0.8), angle=90, size=3) +
	geom_line(data= chr4B_plot_tbl, aes(x = physical_position, y = MAF, color = Improvement_status, group = Accession, alpha= Improvement_status), size=0.3) +
	facet_grid(chromosome~., scales="free") +
	ylim(scan_ylim) +
	xlab(scan_xlab) +
	ylab(scan_ylab) +
	scale_color_manual(values = scan_plot_colours) +
	scale_alpha_manual(values = scan_plot_alpha) +
	scan_plot_theme +
	NULL

################################
################ combine plots
################################

outlier_grob<-ggplotGrob(outlier_plots)
g1 <- ggplotGrob(chr3A_plot + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), strip.background = element_blank(), strip.placement = "outside"))
g2 <- ggplotGrob(chr3B_plot + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), strip.background = element_blank(), strip.placement = "outside"))
g3 <- ggplotGrob(chr4B_plot + theme(legend.position="none", axis.title.x=element_blank(), axis.title.y=element_blank(), strip.background = element_blank(), strip.placement = "outside"))
legend = gtable_filter(ggplot_gtable(
		ggplot_build(chr3A_plot+ guides(
			colour = guide_legend(override.aes = list(alpha = 1, size=1))))), "guide-box")

scans<-grid.arrange(legend, g1, g2, g3 , ncol=1, heights=c(0.3,1,1,1), 
  bottom=textGrob("median physical position of sliding window (mb)", gp=gpar(fontsize=12)),
  left=textGrob(paste0("frequency of minor allele"), gp=gpar(fontsize=12), rot=90))

Alabel<-textGrob("A", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Blabel<-textGrob("B", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))

pdf(file=paste0(plot_export_prefix,"TtBtr1_and_chr4B_all_SNPs.pdf"), height= 4.8*(6/3.3), width=8*1.025)
grid.arrange(Alabel, outlier_grob, Blabel, scans,  nrow=2, widths=c(0.025,1), heights=c(1.5,3.3))
dev.off()


################################
######## Statistics of interest
################################

################################
######## Defining sweep regions
################################

#repeat the dom MAF analysis but use physical window sizes rather than SNP windows:

MAF_by_physical<-MAF_by_physical_window(traw_dom_MAF, window_size=1000000, increments_per_window=2, location)

#Only going to look at domesticated accessions, and only in the vicinity of the genes of interest

#Filtering to the region on chr3A
#the longest continuous window where all domesticated MAF are less than 0.05 is from 59-63mb (recall that rounded gives the midpoint of the sliding window)
max_MAF_by_region(MAF_by_physical, domesticated_columns, "chr3A", lower_bound_mb= chr3A_Br$physical_position-padding_lower, upper_bound_mb = chr3A_Br$physical_position+ padding_upper, threshold=0.05)

#region on chr3B
#The longest *continuous* window is from 94mb-99.5mb
max_MAF_by_region(MAF_by_physical, domesticated_columns, "chr3B", lower_bound_mb= chr3B_Br$physical_position-padding_lower, upper_bound_mb = chr3B_Br$physical_position+ padding_upper, threshold=0.05) %>% print(n=100)

#region on chr4B
#the longest continuous window is from 509-512mb.
max_MAF_by_region(MAF_by_physical, domesticated_columns, "chr4B", lower_bound_mb= 515-50, upper_bound_mb = 515+ 50, threshold=0.05) %>% print(n=100)

################################
######## number of outlier regions within sweeps
################################

#Check the number of UC10164 calls in these regions:
filter(MAF_by_physical, 
	Accession==sample, 
	(chromosome=="chr3A" & rounded/1000000>59 & rounded/1000000<63) |
	(chromosome=="chr3B" & rounded/1000000>94 & rounded/1000000<99.5) |
	(chromosome=="chr4B" & rounded/1000000>509 & rounded/1000000<512)) %>%
	summarise(num_calls_sum=sum(num_calls), minor_alleles_sum=sum(minor_alleles), MAF=minor_alleles_sum/num_calls_sum)

#for comparison, can look at the average number of minor alleles in these regions across all modern domesticated accessions. 
dom_sweeps<-ungroup(MAF_by_physical) %>% group_by(Accession) %>%
filter( 
	Accession %in% domesticated_columns, 
	(chromosome=="chr3A" & rounded/1000000>59 & rounded/1000000<63) |
	(chromosome=="chr3B" & rounded/1000000>94 & rounded/1000000<99.5) |
	(chromosome=="chr4B" & rounded/1000000>509 & rounded/1000000<512)) %>%
	summarise(num_calls_sum=sum(num_calls), minor_alleles_sum=sum(minor_alleles), MAF=minor_alleles_sum/num_calls_sum)
summary(dom_sweeps$MAF)

################################
######## Implied selection coefficients from sweep sizes
################################

#From Olsen et al. (2006). "There are several methods to estimate the strength of selection on a gene, s, from molecular population genetic data. One such approach utilizes the reduced variation surrounding a target of directional selection that arises from genetic hitchhiking of neutral sites linked to the selected site. The distance between a selected and hitchhiking neutral site (d) is dependent on the strength of selection (either as the selection coefficient s or as the population selection parameter α = 2Nes) and on the recombination rate (c or as the population recombination parameter C = 2Nec) in this genomic region. The relationship among these variables is d = 0.01α/C = 0.01s/c (Kaplan et al. 1989),"
#From Kaplan, Hudson, and Langley: "This means that the expected level of variation will be substantially reduced for all the sites within a physical distance of (0.01)alpha/C base pairs of a locus at which a selected substitution has recently occurred. 
#"We can calculate the strength of selection using either direct estimates of recombination or the population recombination parameter. Using crossover data between wx mutants, a direct estimate of the intragenic recombination frequency for Wx has been estimated previously as c = 3.67 × 10−7/bp (Inukai et al. 2000). If we assume the span of the affected genomic region to be 250 kb, then the mean distance between the selected site and the farthest hitchhiking site would be 125 kb, and the inferred selection coefficient would be ∼4.59."
#From Walsh and Lynch (2018):
#"The sweep signature spans 250kb, encompassing ~40 genes. Further, there is a strong EHH signal (Table 9.3) around Waxy, and alleles from temperate japonica lines show a highly negative Tajima's D. Oslon et al. assumed that c=3.7x10^(-7) per bp (Inukai et al. 2000) and used Equation 8.6b to estimate the strength of selection as
#s=3.7*10^(-7)*250000/0.02=4.6
#"This estimated value implies incredibly strong selection, with individuals carrying this alleles leaving (on average) close to five times as many offspring as those without it. However, this estimate does not account for the reduction in recombination from selfing. Using the effective recombination rate (Equation 9.44) reduces the estimate to a more modest value of s ~0.1 (assuming a high selfing rate of eta=0.99)
#eq 9.44: c*~c(1-eta/(2-eta))

#We have map distances in cM. Convert to recombination rates (in crossover units, c.u.) using Haldane's map function:
#Pr(recombination| linkage of x cM)= (1-exp(-2x/100))/2
#we have x in units of cM per mb, so divide by 10^6
xparam=0.2/1000000
cparam=(1-exp(-2*xparam/100))/2
#However, we are going to adjust c to cstar, which accounts for selfing
etaparam=0.99
cstarparam=cparam*(1-etaparam/(2-etaparam))
#if d=0.01s/c, s=c*d/0.02 (0.02 accounts for both directions of the sweep)
#distance of 4mb
dparam=4*1000000
sel=cstarparam*dparam/0.02

calculate_sel<-function(xparam_in_cM_per_mb, eta_param, dparam){
	cparam=(1-exp(-2* xparam_in_cM_per_mb/(100*10^6)))/2 #recombination rate from cM genetic map distance
	cstarparam=cparam*(1-eta_param/(2-eta_param)) #effective recombination rate given selfing rate
	sel=cstarparam*dparam/0.02 #selection parameter
	return(sel)
}

calculate_sel(0.1, 0.99, 4*1000000)
calculate_sel(0.1, 0.995, 4*1000000)
calculate_sel(0.5, 0.99, 4*1000000)
calculate_sel(0.5, 0.995, 4*1000000)
calculate_sel(0.1, 0.99, 5.5*1000000)
calculate_sel(0.1, 0.995, 5.5*1000000)
calculate_sel(0.5, 0.99, 5.5*1000000)
calculate_sel(0.5, 0.995, 5.5*1000000)

calculate_sel(0.05, 0.99, 3*1000000)
calculate_sel(0.05, 0.995, 3*1000000)
calculate_sel(0.2, 0.99, 3*1000000)
calculate_sel(0.2, 0.995, 3*1000000)

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

}