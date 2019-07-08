library(ggplot2)
library(gridExtra)
library(gridGraphics)

root_tree<-function(tree){
	rooted<-root(tree, outgroup = "Outgroup", resolve.root = TRUE)
	#Set the root position halfway along the outgroup edge
	root_edge<-length(rooted$edge.length)
	len<-rooted$edge.length[root_edge]/2
	rooted$edge.length[1]<-len
	rooted$edge.length[root_edge]<-(rooted$edge.length[root_edge]-len)
	return(rooted)
}

adjust_node_labels<-function(tree){
	#remove node labels where bootstrap support=100 (and root name)
	tree$node.label[tree$node.label=="Root"]<-NA
	tree$node.label[as.numeric(tree$node.label)==100]<-NA
	return(tree)
}