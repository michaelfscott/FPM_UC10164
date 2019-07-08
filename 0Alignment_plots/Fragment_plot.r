source("Fragment_plot_functions.r")

#parameters
S1_path="../../4Alignments/1Align_Emmer_Collapsed/Merged/UC10164S1/"
S2_path="../../4Alignments/1Align_Emmer_Collapsed/Merged/UC10164S2/"
S1_path_bread="../../4Alignments/1Align_Bread_Collapsed/Merged/UC10164S1/"
S2_path_bread="../../4Alignments/1Align_Bread_Collapsed/Merged/UC10164S2/"
min_mapQs<-c(0,1,10,20,25,30,35)
columnNames<-c("chromosome", "start", "end", "raw", "dedup")

#get data 
lengths<-get_fragment_lengths(S1_path, S2_path)
mis<-get_misincorporation(S1_path, S2_path) %>%
  filter(Position>=(-20) & Position <=20)
tbl_bread<-get_all_coverage_chromosomal(S1_path_bread, S2_path_bread, columnNames, min_mapQs, by_subgenome=TRUE) %>%
  filter(status=="raw", !subgenome=="Un")
tbl_bread_mapped<-get_perc_coverage_chromosomal(S1_path_bread, S2_path_bread, columnNames, min_mapQs) %>%
  filter(status=="raw") 
tbl_emmer<-bind_rows(list(
	SNP_sites= get_all_coverage_SNP(S1_path, S2_path, columnNames, min_mapQs) , 
	all= get_all_coverage_chromosomal(S1_path, S2_path, columnNames, min_mapQs)), .id="location") %>%
	filter(status=="raw")

library(gridGraphics)
library(gridExtra)

length_grob<-ggplotGrob(plot_lengths(lengths))
mis_grob<-ggplotGrob(plot_mis(mis))
bread_cov_grob<-ggplotGrob(plot_bread_AB_cov(tbl_bread))
bread_perc_grob<-ggplotGrob(plot_bread_subgenome_perc(tbl_bread_mapped))
emmer_grob<-ggplotGrob(plot_emmer_full_SNP(tbl_emmer))

Alabel<-textGrob("A", x = unit(0, "npc"), y= unit(1, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Blabel<-textGrob("B", x = unit(0, "npc"), y= unit(1, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Clabel<-textGrob("C", x = unit(0, "npc"), y= unit(1, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))
Dlabel<-textGrob("D", x = unit(0, "npc"), y= unit(1, "npc"), just=c("left","top"), gp=gpar(col="black", fontsize=18, fontface="bold"))

hlay <- rbind(c(1,2,2,3,4),
              c(1,1,1,3,3),
              c(5,6,7,8,9))

plot_export_prefix="plots/"
pdf(file=paste0(plot_export_prefix, "Fragment_alignments.pdf"), height=9, width=(9))
	grid.arrange(Alabel, mis_grob, Blabel, length_grob, Clabel, bread_cov_grob, bread_perc_grob, Dlabel, emmer_grob, heights=c(1,0.1,1.2), widths=c(0.1,0.7,0.8,0.1,1), layout_matrix=hlay)
dev.off()

