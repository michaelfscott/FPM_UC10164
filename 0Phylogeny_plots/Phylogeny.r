args <- commandArgs(trailingOnly = TRUE)
if (!length(args)==3) {
  stop("supply arguments: (1) metadata_filename (2) tree_filename (3) export_dir")
} else {

#metadata_filename="../../2Reference_Genomes/Zavitan_v2/metadata/origins.txt"
metadata_filename=args[1]
#tree_filename="../../6SNPhylo/2SNPhylo_result_bedfilter_no_indels_altref_altbam/mapQ30/excludeS1/snphylo.output.bs.tree"
tree_filename=args[2]
#plot_export_prefix="plots/"
plot_export_prefix=args[3]

source("../0Rfunctions/functions.r")
source("Phylogeny_plot_functions.r")

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
rotated <-adjust_node_labels(rooted) #remove node labels where bootstrap support=100 (and root name)

quartz(height=9, width=5)
plot.phylo(rotated, 
	show.node.label=FALSE, 
	tip.col= sapply(as.numeric(as.factor(location[,"Region"])), function(x)(customPalette[-1][x])), 
	cex=0.75, edge.width =2, root.edge=TRUE, label.offset=0.01, 
	rotate.tree=90, no.margin=TRUE, 
	adj=1)

#add bootstrap support
node_numbers<-as.numeric(which(!is.na(rotated$node.label)))+length(tree$tip.label) #have to add to take account of tip nodes
node_text<-rotated$node.label[!is.na(rotated$node.label)]
nodelabels(text= node_text, node=node_numbers,cex=0.5, bg="white", frame="circle")

phylo_grob<-grab_grob()
dev.off()

pdf(file=paste0(plot_export_prefix, "Tree_excludeS1.pdf"), height=9, width=(5))
grid.arrange(phylo_grob)
dev.off()

}
