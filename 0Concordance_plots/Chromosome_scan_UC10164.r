args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==3) {
  stop("supply arguments: (1) extra_data_dir (2) traw_file (3) plot_export_prefix")
} else {

source("../0Rfunctions/functions.r")
source("../0Rfunctions/dp.r")


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

traw<-read_tsv(traw_file) %>% join_chromosome_halves_rename_samples() %>% select(-exclude_samples)
location<-load_metadata(metadata_filename) %>% 
	filter(Accession %in% colnames(traw) & !Accession %in% exclude_samples) 

###########
# similarity plots for UC10164
###########

snp_window_size<-50
window_move_size<-25

sample="UC10164S1"

sim<-similarity_by_SNP_window(traw, sample, snp_window_size, window_move_size, location)

left_join(sim, location, by="Accession") %>%
ggplot(aes(colour=Region, x = physical_position, y = MAF)) + 
  geom_line(aes(group = Accession), size=0.1) +
  ylim(-0.1,1) +
  facet_grid(chromosome~.) +
  theme(legend.position="top") +
  xlab("median physical position of window (mb)") +
  ylab(paste0("differences from ", sample, " (in ", snp_window_size, " marker sliding window)")) +
#  geom_point(data=mosaic_df, aes(x=physical_position, y=-0.1, colour= haplotype_region) , pch="|") + 
  scale_color_manual(values = customPalette) +
  facet_grid(chromosome~.) +
  NULL

###########
# Dynamic Programming 
###########

#don't allow het calls or NA calls in sample and remove sample from test
raw<-t(traw[,7:ncol(traw)])
raw<-raw[location$Accession,]
test<-raw[,! raw[sample,]==1 & !is.na(raw[sample, ])]
#convert all scores such that the sample of interest has value 0
test<-apply(test,2,function(x)(if(x[sample]==0){x}else{abs(x-2)}))
#remove sample from test
test<-test[-which(rownames(test)==sample),]

map_to_use<-data.frame(traw[!traw[,sample]==1 & ! is.na(traw[, sample]),c("chromosome", "position")])

source("dp.R")

jumpcost=2.5
threshold=0.075

mosaic<-haplotype.mosaic(test, map_to_use, jumpcost=jumpcost, threshold=threshold)

mosaic_chr<-rep(names(mosaic), sapply(mosaic, nrow))
mosaic_pos<-do.call(rbind,mosaic)[,2]/1000000 #convert to mb
mosaic_haplotype<-do.call(rbind,mosaic)[,3]
#Re-order so that accession match between genotypes and metadata 
sorted_regions<-left_join(data.frame(Accession=rownames(raw)), location)$Region
mosaic_haplotype_region<-sorted_regions[mosaic_haplotype]
mosaic_haplotype_region<-gsub("Ancient_Egyptian", " No Good Match", mosaic_haplotype_region)
mosaic_df<-data.frame(chromosome= mosaic_chr, physical_position= mosaic_pos, haplotype=mosaic_haplotype, haplotype_region=mosaic_haplotype_region)

pdf(file=paste0(plot_export_prefix,"/chromosome_scan_JC",jumpcost,"_TH",threshold,".pdf"), height=15, width=8)
filter(sim, !Accession == sample) %>%
left_join(location, by="Accession") %>%
ggplot(aes(colour=Region, x = physical_position, y = MAF)) + 
  geom_line(aes(color = Region, group = Accession), size=0.1) +
  ylim(-0.1,1) +
  facet_grid(chromosome~.) +
  theme(legend.position="top") +
  xlab("median physical position of window (mb)") +
  ylab(paste0("differences from ", sample, " (in ", snp_window_size, " marker sliding window)")) +
  geom_point(data=mosaic_df, aes(x=physical_position, y=-0.1, colour= haplotype_region) , pch="|") + 
  scale_color_manual(values = customPalette) +
  facet_grid(chromosome~.) +
  NULL
dev.off()

