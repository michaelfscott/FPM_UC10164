library(ggplot2)
library(ggfortify)
library(gplots)
library(stringr)
library(reshape2)
library(ape)
library(dplyr)
library(readr)

#not sure how many are necessary
library(ggmap)
library(maps)
library(mapdata)
library(grid)
library(gtable)
library(gridExtra)

plot_export_prefix="plots/"
input_prefix="../../5Variants/3Modern_Samples_Merge_bedfilter_noindels_altref_altbam/mapQ30/S1minDP2/merged"

traw<-read_tsv(paste0(input_prefix,".traw"))
raw<-t(traw[,7:ncol(traw)])
name=rownames(raw)
rownames(raw)<-ifelse(grepl("Mt_", name)| grepl("Iraq_", name), paste0(word(name,1,sep="_"),"_", word(name,2,sep="_")), word(name,1,sep="_"))
#remove Outgroup for this analysis
raw<-raw[-which(rownames(raw)=="Outgroup"),]

################################
################ rearrange & parameters
################################

location<-read.table("../../2Reference_Genomes/Zavitan_v2/metadata/origins.txt", stringsAsFactors=FALSE, header=TRUE) %>% 
	rbind(c("UC10164S1","Egypt","UC10164S1", "dicoccum", "domesticated", "Ancient_Egyptian", 26.93, 31.48)) %>%
	rbind(c("Outgroup", NA, "Outgroup", "Outgroup", "Outgroup", "Outgroup", NA,NA))
location$East<-as.numeric(location$East)
location$North<-as.numeric(location$North)
rownames(location)<-location$Accession
location<-location[rownames(location) %in% rownames(raw),]
location<-location[rownames(raw),]
#relabel "Region" of two samples that seem to fall into the wrong group genetically and geographically
location[location$Label %in% c("WE-19", "WE-20"),"Region"]<-"Southern_Levant"

#create "type" of both improvment status and region
location$type=paste(location$Improvement_status, location$Region, sep="_")

#cbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
#customPalette <- c("#000000", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#0072B2", "#56B4E9")
customPalette <- c("#000000", "#E69F00", "#D55E00", "#009E73", "#CC79A7", "#000066", "#56B4E9")
label_names<- c("Ancient Egyptian", "Caucasus (Dom)", "Eastern Europe (Dom)", "Indian Ocean (Dom)", "Mediterranean (Dom)", "Northern Levant (Wild)", "Southern Levant (Wild)")

replace_na_with_median<-function(vector){
vector[is.na(vector)]<-median(vector, na.rm=TRUE)
return(vector)
}

################################
################ maps
################################

world <- map_data("world") 

longlat1<-as.numeric(c(min(location$East)-5, max(location$East)+5, min(location$North)-5, max(location$North)+5))
longlat2<-c(34.27, 44.75, 31.20, 39.38)
#as.numeric(c(min(location$East[location$Region=="Northern_Levant"])-1, max(location$East[location$Region=="Northern_Levant"])+1, min(location$North[location$Region=="Northern_Levant"])-1, max(location$North[location$Region=="Northern_Levant"])+1.5))
pol<-data.frame(xmin= longlat2[1]-0.5,xmax=longlat2[2]+0.5 ,ymin=longlat2[3]-0.375 ,ymax=longlat2[4]+0.375)

labs<-data.frame(names=location$Label, long=location$East, lat=location$North, type=location$type)

zoomout <-ggplot() + geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey70", color="grey70")+ 
  geom_point(data = labs, aes(x = long, y = lat, shape=type, colour= type), size=1.25, stroke=1.25, alpha=0.75) +
#  geom_text(data = labs, aes(x = long, y = lat, colour= factor(region), label=names), hjust=1, vjust=0) +
#  scale_color_manual(values= cbPalette[-4]) +
  scale_shape_manual(labels= label_names , values=c(19,19,19,19,19,21, 21))+
  scale_color_manual(labels= label_names , values= customPalette) +
  coord_fixed(xlim = longlat1[1:2],  ylim = longlat1[3:4], ratio = 1.3) +
  geom_rect(data = pol, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha=0, colour="black", size = 0.75, linetype=1) +
  scale_y_continuous(expand = c(0,0)) + scale_x_continuous(expand = c(0,0)) +
  theme_nothing() +
#  ggsn::scalebar(data = world, dist=100, location="topright", st.size=2) +
#  theme(axis.text.x =element_blank(),axis.text.y= element_blank(), axis.ticks=element_blank(),axis.title.x =element_blank(),axis.title.y= element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none")
  NULL

zoomin<-ggplot() + 
  geom_rect(data = pol, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), stroke=2, fill="white", size = 0, linetype=1) +
  geom_polygon(data = world, aes(x=long, y = lat, group = group), fill="grey70", color="grey70")+ 
  geom_point(data = labs, aes(x = long, y = lat, shape=type, colour= type), size=1.25, stroke=1.25, alpha=0.75) +
#  geom_text(data = labs, aes(x = long, y = lat, colour= factor(region), label=names), hjust=1, vjust=0) +
  scale_shape_manual(labels= label_names , values=c(19,19,19,19,19,1, 1))+
  scale_color_manual(labels= label_names , values= customPalette) +
  coord_fixed(xlim = longlat2[1:2],  ylim = longlat2[3:4], ratio = 1.3) +
  geom_rect(data = pol, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), alpha=0, colour="black", size = 1.5, linetype=1) +
  theme_nothing()
#  theme(axis.text.x =element_blank(),axis.text.y= element_blank(), axis.ticks=element_blank(),axis.title.x =element_blank(), axis.title.y= element_blank(), panel.background = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position="none", panel.border = element_blank(), panel.spacing=element_blank())

v1<-viewport(width = 1, height = 1, x = 0.5, y = 0.5) #plot area for the main map
v2<-viewport(width = 0.5/((longlat2[1]-longlat2[2])/(longlat2[3]-longlat2[4])), height = 0.5, x = 1-(0.275)/2, y = 1-(0.5)/2) #plot area for the inset map

#print(zoomout,vp=v1) 
#pdf(file=paste0(map_export_prefix,"map_with_zoom.pdf"), height=4, width=3.1*(longlat1[1]-longlat1[2])/(longlat1[3]-longlat1[4]))
#zoomout
#plot(zoomin,vp=v2, add=TRUE)
#dev.off()

################################
################ PCA non-ancient samples and then project
################################

#PCA, first make replace all NA's with the median. 
sample="UC10164S1"
modern_raw<-raw[-which(rownames(raw)==sample),]
test_modern<-apply(modern_raw,2,replace_na_with_median)
prin_comp_modern<-prcomp(test_modern, scale.=TRUE)
#predict this sample using the prin_comp of the modern ones
predict_prin_comp<-predict(prin_comp_modern, t(data.frame(raw[which(rownames(raw)==sample),])))

#warning, this step only works because sample is the first one in location
PCAobject_modern=data.frame(location, rbind(predict_prin_comp, prin_comp_modern$x))
eigs <- prin_comp_modern$sdev^2
perc_var_modern<-round(eigs / sum(eigs) *100, digits=2)

#PCA, first make replace all NA's with the median. 
#get PCA for domesticated samples:
modern_raw_dom<-raw[!rownames(raw)==sample & location$Improvement_status=="domesticated",]
test_modern_dom<-apply(modern_raw_dom,2,replace_na_with_median)
test_modern_dom2<-test_modern_dom[,!apply(test_modern_dom,2,function(x)(var(x)==0))]
prin_comp_modern_dom<-prcomp(test_modern_dom2, scale.=TRUE)
#predict this sample using the prin_comp of the modern ones
predict_prin_comp_modern_dom<-predict(prin_comp_modern_dom, t(data.frame(raw[which(rownames(raw)==sample),!apply(test_modern_dom,2,function(x)(var(x)==0))])))

#warning, this step only works because sample is the first one in location
PCAobject_modern_dom=data.frame(location[location$Improvement_status=="domesticated",], rbind(predict_prin_comp_modern_dom, prin_comp_modern_dom$x))
eigs <- prin_comp_modern_dom$sdev^2
perc_var_modern_dom<-round(eigs / sum(eigs) *100, digits=2)

################################
################ PCA plots
################################

d1<-data.frame(PCAobject_modern[,1:10], y= PCAobject_modern[,"PC2"], panel="PC2", samples="all")
d2<-data.frame(PCAobject_modern[,1:10], y= PCAobject_modern[,"PC3"], panel="PC3", samples="all")
d1_dom<-data.frame(PCAobject_modern_dom[,1:10], y= PCAobject_modern_dom[,"PC2"], panel="PC2", samples="domesticated")
d2_dom<-data.frame(PCAobject_modern_dom[,1:10], y= PCAobject_modern_dom[,"PC3"], panel="PC3", samples="domesticated")
d<-rbind(d1 , d2)

p<-ggplot(data=d, scale.y="free", aes(x=PC1,y=y, shape=type, col=type)) +
  geom_point(size=1.25, stroke=1.25, alpha=0.75)+ 
  scale_shape_manual(labels= label_names , values=c(19,19,19,19,19,1, 1))+
  scale_color_manual(labels= label_names , values = customPalette) +
  xlab(paste0("PC1 (", perc_var_modern[1], "%)"))+
  theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) +
#  coord_fixed() +
  facet_grid(panel~ samples,scales="free", switch="both", 
  	labeller = as_labeller(c(all = paste0("PC1 (", perc_var_modern[1], "%)") , PC2 = paste0("PC2 (", perc_var_modern[2], "%)"), PC3 = paste0("PC3 (", perc_var_modern[3], "%)")) )) +
  theme(strip.background = element_blank(),
	strip.placement = "outside") +
  NULL

d_dom<-rbind(d1_dom , d2_dom)
p_dom<-ggplot(data= d_dom, scale.y="free", aes(x=PC1,y=y, shape= type, col= type)) +
  geom_point(size=1.25, stroke=1.25, alpha=0.75)+ 
  scale_shape_manual(labels= label_names , values=c(19,19,19,19,19,1, 1))+
  scale_color_manual(labels= label_names , values = customPalette) +
  xlab(paste0("PC1 (", perc_var_modern[1], "%)"))+
  theme(legend.title=element_blank(), axis.title.x=element_blank(), axis.title.y=element_blank()) +
  facet_grid(panel~ samples,scales="free", switch="both", 
  	labeller = as_labeller(c(domesticated = paste0("PC1 (", perc_var_modern_dom[1], "%)") , PC2 = paste0("PC2 (", perc_var_modern_dom[2], "%)"), PC3 = paste0("PC3 (", perc_var_modern_dom[3], "%)")) )) +
  theme(strip.background = element_blank(),
	strip.placement = "outside") +
#  coord_fixed() +
  NULL

################################
################ Unique windows plot
################################

threshold=0.05
load(file= paste0("../../5Compare/scripts/RData/conc_",threshold,".RData"))
load(file= paste0("../../5Compare/scripts/RData/ncol_res.RData"))
levels(as.factor(tmp$Type))
tmp<-data.frame(location, all_low_conc= all_low_conc) %>% tbl_df() %>%
  mutate(Type=paste(Improvement_status, Region, sep="_"))
unique_plot<-ggplot(tmp, aes(y= all_low_conc/ncol_res, x=factor(tmp$Region, levels=c("Ancient_Egyptian", "Indian_Ocean", "Mediterranean", "Eastern_Europe", "Caucasus", "Northern_Levant", "Southern_Levant")), colour=Type, pch=Type)) +
  geom_jitter(width=0.25, size=1.25, stroke=1.25, alpha=0.75) +
  scale_colour_manual(labels= label_names , values = customPalette) +
#  scale_alpha_discrete( range=c(0.75,0.75,0.75,0.75,1,1)) +
  scale_shape_manual(labels= label_names , values=c(19,19,19,19,19,1, 1))+
  theme(legend.position="none", 
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank()) + 
  labs(y=paste0("Fraction \'unique\' haplotypes")) +
  NULL

################################
################ Combine
################################

g_unique<-ggplotGrob(unique_plot)
g <- ggplotGrob(p + theme(legend.position="none"))
g_dom <- ggplotGrob(p_dom + theme(legend.position="none"))
#get the legend as a separate grob
legend = gtable_filter(ggplot_gtable(ggplot_build(p+ theme(legend.text=element_text(size=10) ) + guides(
			colour = guide_legend(nrow=2), pch= guide_legend(nrow=2)))), "guide-box")

#convert zoomout to a grob
g_map<-ggplotGrob(zoomout)

#plot area for the inset map, which will be plotted onto the whole figure:
v2<-viewport(width = 0.35/((longlat2[1]-longlat2[2])/(longlat2[3]-longlat2[4])), height = 0.35, x = 1-(0.285)/2, y = 1-(0.48)/2) 

Alabel<-textGrob("A", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Blabel<-textGrob("B", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Clabel<-textGrob("C", x = unit(0, "npc"), y= unit(0.95, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))

hlay <- rbind(c(1,1,1,1,1),
              c(2,3,3,3,3),
              c(4,5,6,7,8))

pdf(file=paste0(plot_export_prefix, "Map_and_PCA_and_unique.pdf"), height=7*1.2, width=7*1.05)
grid.arrange(legend, Alabel, g_map,Blabel, g, g_dom, Clabel, g_unique, heights=c(0.2,1,1), widths=c(0.025,0.33,0.33,0.025,0.34), layout_matrix=hlay)
plot(zoomin,vp=v2, add=TRUE)
dev.off()

pdf(file=paste0(plot_export_prefix, "Map_and_PCA_and_unique_tmp.pdf"), height=7*1.2, width=7*1.05)
grid.arrange(legend, Alabel, g_map,Blabel, g, g_dom, Clabel, Clabel, heights=c(0.2,1,1), widths=c(0.025,0.33,0.33,0.025,0.34), layout_matrix=hlay)
plot(zoomin,vp=v2, add=TRUE)
dev.off()
