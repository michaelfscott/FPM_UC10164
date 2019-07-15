args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==4) {
  stop("supply arguments: (1) metadata_filename (2) tree_filename (3) Dstats_path (4) export_dir")
} else {

#metadata_filename="../../2Reference_Genomes/Zavitan_v2/metadata/origins.txt"
metadata_filename=args[1]
#tree_filename="../../6SNPhylo/2SNPhylo_result_bedfilter_no_indels_altref_altbam/mapQ30/S1minDP2/incl_UC10164S1/snphylo.output.bs.tree"
tree_filename=args[2]
#Dstats_path="../../8AdmixTools/2AdmixTools_out_bedfilter_no_indels_altref_altbam/mapQ30/S1minDP2/"
Dstats_path=args[3]
#plot_export_prefix="plots/"
plot_export_prefix=args[4]

source("../0Rfunctions/functions.r")
source("../0Rfunctions/Phylogeny_plot_functions.r")

######
#read data
######

location<-load_metadata(metadata_filename)
tree<-read.tree(tree_filename)
#adjust location to have only the appropriate entries and re-ordered
location<-left_join(data.frame(Accession=tree$tip.label), location)

customPalette <- c("#000000", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#000066", "#000000", "#56B4E9")

######
#plot tree
######

rooted<-root_tree(tree)
#rotate to show introgression between the "Tabhinga group" and UC10164
rotated<-rotate(rotate(rooted, 69), 110) %>% adjust_node_labels() #remove node labels where bootstrap support=100 (and root name)

quartz(height=7, width=7)
plot.phylo(rotated, 
	type="fan", show.node.label=FALSE, 
	tip.col= sapply(as.numeric(as.factor(location[,"Region"])), function(x)(customPalette[x])), 
	cex=0.75, edge.width =2, root.edge=TRUE, label.offset=0.01, 
	rotate.tree=90, no.margin=TRUE, 
	adj=1, y.lim=c(-0.5,0.4), x.lim=c(-0.8,0.8)) 

#add bootstrap support
node_numbers<-as.numeric(which(!is.na(rotated$node.label)))+length(tree$tip.label) #have to add to take account of tip nodes
node_text<-rotated$node.label[!is.na(rotated$node.label)]
nodelabels(text= node_text, node=node_numbers,cex=0.5, bg="white", frame="circle")

phylo_grob<-grab_grob()
dev.off()

########
###Dstat results
########

column_names<-c("result", "Outgroup","SL_Accession","IO","UC10164", "D_stat", "SE", "Z", "BABA","ABBA","num_SNPs")
Dstat<-read.table(paste0(Dstats_path,"Split_Southern_Levant_qpDstat.5mb.results.txt"), stringsAsFactors=FALSE, skip=40, col.names=column_names)

rownames(Dstat)<-Dstat$SL_Accession
SL_ordering<-c( "PI487264", "PI428132", "Mt_Hermon", "PI352322", "PI466957", "PI538673", "PI471021", "Gamla", "Qazerin", "PI466946", "PI428122", "PI471038", "Mt_Gerizim", "Nesher", "PI471060","Tabigha", "PI467008", "PI466950", "PI470962")

Dstat_plot<-ggplot(Dstat, aes(x= D_stat, y= factor(SL_Accession, levels= SL_ordering))) + 
    geom_errorbarh(aes(xmin= D_stat-2*SE, xmax= D_stat+2*SE), height=.1) +
    geom_point() +
    geom_vline(xintercept=0, linetype="dashed") +
    xlab("D statistic (+/- 2SE)") +
    ylab("Wild Southern Levant Accession") +
    theme_minimal() + theme(axis.text.y = element_text(colour = customPalette[8]),
		panel.grid.minor = element_blank())

Dstat_Grob <- ggplotGrob(Dstat_plot)

########
### example tree
########

tree <- read.tree(text = "(Outgroup,(Southern_Levant,(Indian_Ocean,UC10164)));")
quartz(height=7*0.25, width=7*0.2)
plot.phylo(tree, type="cladogram", direction="upwards", tip.col=customPalette[c(1,8,4,7)], label.offset=0.1, edge.width =2, cex=0.8, no.margin=TRUE, x.lim=c(-1,6))
example_grob<-grab_grob()
dev.off()

########
### combine plots
########

Alabel<-textGrob("A", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Blabel<-textGrob("B", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))

hlay <- rbind(c(5,1,6,NA, 2,3),
              c(5,1,6,4,4,4))

pdf(file=paste0(plot_export_prefix, "Tree_and_Dstat_with_example_tree.pdf"), height=7, width=(7*1.5))
grid.arrange(phylo_grob, example_grob, example_grob, Dstat_Grob, Alabel, Blabel, heights=c(0.25,0.7), widths=c(0,1,0,0.1,0.2,0.2), layout_matrix=hlay)
dev.off()

########
###all Dstat results table
########

Dstat_100kb<-read.table(paste0(Dstats_path,"/Split_Southern_Levant_qpDstat.100kb.results.txt"), stringsAsFactors=FALSE, skip=40, col.names=column_names) %>%
transmute(SL_Accession=SL_Accession, SE_100kb=SE, Z_100kb=Z)
Dstat_500kb<-read.table(paste0(Dstats_path,"/Split_Southern_Levant_qpDstat.500kb.results.txt"), stringsAsFactors=FALSE, skip=40, col.names=column_names) %>%
transmute(SL_Accession=SL_Accession, SE_500kb=SE, Z_500kb=Z)
Dstat_1mb<-read.table(paste0(Dstats_path,"/Split_Southern_Levant_qpDstat.1mb.results.txt"), stringsAsFactors=FALSE, skip=40, col.names=column_names) %>%
transmute(SL_Accession=SL_Accession, SE_1mb=SE, Z_1mb=Z)
Dstat_10mb<-read.table(paste0(Dstats_path,"/Split_Southern_Levant_qpDstat.10mb.results.txt"), stringsAsFactors=FALSE, skip=40, col.names=column_names) %>%
transmute(SL_Accession=SL_Accession, SE_10mb=SE, Z_10mb=Z)
Dstat_50mb<-read.table(paste0(Dstats_path,"/Split_Southern_Levant_qpDstat.50mb.results.txt"), stringsAsFactors=FALSE, skip=40, col.names=column_names) %>%
transmute(SL_Accession=SL_Accession, SE_50mb=SE, Z_50mb=Z)

output<-full_join(Dstat, Dstat_100kb) %>%
full_join(Dstat_500kb) %>%
full_join(Dstat_1mb) %>%
full_join(Dstat_10mb) %>%
full_join(Dstat_50mb) %>%
mutate(Z_5mb=Z, SE_5mb=SE) %>%
select("SL_Accession", "BABA", "ABBA", num_SNPs, D_stat, SE_100kb, SE_500kb, SE_1mb, SE_5mb, SE_10mb, SE_50mb, Z_100kb, Z_500kb, Z_1mb, Z_5mb, Z_10mb, Z_50mb) %>%
arrange(factor(SL_Accession, levels= rev(SL_ordering)))

write.csv(output, file= paste0(plot_export_prefix, "Dstats_table.csv"), row.names=FALSE, quote=FALSE)

}
