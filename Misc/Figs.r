# Information  --------------------------------------------------------------------
## Script name: Fig.R
## Description: Figures scratch script.
## 
# Body  --------------------------------------------------------------------
source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
RVMT <- "/media/HDD1/uri/RNA_Vir_MTs/RVMT/"
RNA_Vir_MTs_dir <- "/media/HDD1/uri/RNA_Vir_MTs/V3/"
versiza <- "Wolf"
current_path <- p0(RNA_Vir_MTs_dir, versiza, "/")
setwd(current_path)
library(ape)
library(castor)
library(phangorn)

THREADS <- 11
Memory <- 11800

#### Read Info ####
IDFT <- Rename1Col(fread0("../Vfin/Metadata/IDFT.14092021.tsv"), "Contig", "ND")
IDFT.fasta <- readDNAStringSet("/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ALL_nuc_3007.fasta")
WorkTree <- as.phylo(read.tree("./210822/ali.210822m/210822m.77510.09.lab_mod.tre"))
AllRdRps = readAAStringSet("../Vre_3/Double_Clustering3/Hybrid_RdRps.faa")

# Figures playground  --------------------------------------------------------------------
library(readr)
library(Biostrings)
library(dplyr)
library(plotly)
library(reshape2)
library(ggplot2)
library(hrbrthemes)
library(viridis)


##### RCR90/RvANI90 Cluster size distribution ####
DFFigS4a <- plyr::count(plyr::count(SimonDFT2$RCR90[wana(SimonDFT2$RCR90)])$freq)
FigS4a <- ggplot(DFFigS4a) + 
  geom_histogram(aes(x=freq),color="black",bins = 500,) + scale_x_log10()+ scale_y_log10()

FigS4a <- ggplot(DFFigS4a) + 
  geom_point(aes(x=x, y=freq),color="black",alpha=0.5,bins = 500,) + scale_x_log10()+ scale_y_log10()+
  geom_line(aes(x=x, y=freq)) +
  xlab("#Of clusters") + ylab("#Cluster size")  #+ ggtitle(label ="RCR90 Cluster size distribution")


DFFigS4b <- plyr::count(plyr::count(SimonDFT2$RvANI90[wana(SimonDFT2$RvANI90)])$freq)
FigS4b <- ggplot(DFFigS4b) + 
  geom_histogram(aes(x=freq),color="black",bins = 500,) + scale_x_log10()+ scale_y_log10()

FigS4b <- ggplot(DFFigS4b) + 
  geom_point(aes(x=x, y=freq),color="black",alpha=0.5,bins = 500,) + scale_x_log10()+ scale_y_log10()+
  geom_line(aes(x=x, y=freq)) +
  xlab("#Of clusters") + ylab("#Cluster size")  


FigS4 <- egg::ggarrange(plots = list("RCR90" = FigS4a,"RvANI90" = FigS4b),nrow=1,labels = c("RCR90","RvANI90"))# +
  # ggtitle(label ="Cluster size distribution")

pdf(file = "FigS1.pdf",title = "Fig. S1", width = 7.5,height = 5)  
FigS4
dev.off()

##### Domain matches distribution ####
DFFigS5 <- (plyr::count(tmpAllDFs2$New_Name.y[wana(tmpAllDFs2$New_Name.y)]))
FigS5 <- ggplot(DFFigS5) + 
  geom_bar(aes(y=freq,x=x),color="black",stat = "identity") +  scale_y_log10() +
  xlab("#Of clusters") + ylab("#Cluster size")  + ggtitle(label ="Figure S4 - Cluster size distribution")

FigS5 <- ggplot(DFFigS5, aes(y = reorder(x, freq), x = freq)) +
  geom_bar(stat = "identity") +   scale_x_log10() +
  theme(plot.title = element_text(size = 4, hjust = 0.5,face = "bold"))+
  # theme(legend.key.size = unit(x = 0.5, units = "cm"))+
  theme(axis.text.y = element_text(angle = 22,size =1.8))+
  xlab("#Number of HMM matches") + ylab("#Domain function")  + ggtitle(label ="Figure S5 - Domain distribution")


pdf(file = "FigS5.pdf",title = "Fig. S5", width = 4,height = 7.5)  
FigS5
dev.off()

##### Venn diagrams (for novelty visulation) #####
tmpIDFT <- IDFT
Known_Novel_venn <- distinct(tmpIDFT[wana(tmpIDFT$RvANI90),] %>% group_by(RvANI90) %>% summarise(
  HasNovel =  any(Novel),
  HasKnown = any(!Novel)))
VenLst <- list("Novel" = Known_Novel_venn$RvANI90[Known_Novel_venn$HasNovel], "known" = Known_Novel_venn$RvANI90[Known_Novel_venn$HasKnown])
KvsN_ANI90 <- ggvenn::ggvenn(data =  (VenLst))

Known_Novel_venn <- distinct(tmpIDFT[wana(tmpIDFT$RCR90), c("Novel","RCR90")] %>% group_by(RCR90) %>% summarise(
  HasNovel =  any(Novel),
  HasKnown = any(!Novel)))
VenLst2 <- list("Novel" = Known_Novel_venn$RCR90[Known_Novel_venn$HasNovel], "known" = Known_Novel_venn$RCR90[Known_Novel_venn$HasKnown])
KvsN_RCR90 <- ggvenn::ggvenn(data =  (VenLst2))

p123 <- egg::ggarrange(plots = list("Tree Leaves" = KvsN_RCR90,"RvANI90 Clusters" = KvsN_ANI90),nrow=1,labels = c("Tree Leaves","RvANI90 Clusters"))
pdf(file = "Known_vs_Novel_Venn.pdf",title = "Proportion of known, novel, and mixed clusters")
p123
dev.off()


pdf(file = "Known_vs_Novel_pies.pdf",title = "Proportion of known, novel, and mixed clusters")
pie(x = c("Novel" = 65003, "Knwon" = 10540, "Mixed" = 1968))
pie(x = c("Novel" = 117276, "Knwon" = 11537, "Mixed" = 1125))
dev.off()



#### MegaTree but removing the labels ####

p13 <- p12
p13$layers[[14]] <- NULL
QuickRDSSave(object = p13,n.cores = 12,fname = "p13")





