# Information  --------------------------------------------------------------------
## Script name: Cont2leaf.R
## Description: Affiliate contigs (ND) with tree leaves (RCR90).
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

# Body ---------------------------------------------------------------
IDFT_bc = IDFT # Backup.
IDDT <- NovoDFT # Working table, in case "IDFT" carries left overs from previous Cont2leaf attempts.
IDDT$AfLvl <- 42 # Place holder for non-affiliated contigs.
IDDT$node <- NA # Place holder for megatree-node index.
IDDT$RCR90 <- NA # Place holder for megatree-tip label. 

#### Level A. - Megatree leaves (RCR90 representatives) ####
W1 <- wh(IDDT$RID %in% WorkTree$tip.label)
IDDT$AfLvl[W1] = "Lvl A. - Megatree leaves."
IDDT$node[W1] <-  unlist(NLabel2NID(WorkTree,IDDT$RID[W1]))
IDDT$RCR90[W1] <- IDDT$RID[W1]

#### Level B. - best BLASTP hit from RdRp to {0}; pident >=90% AND qcov >=75% ####
BlastpOut = fread("../Wolf/tmp.07102021.50.50.tab",nThread = 10 ,header = F,col.names = c("qid", "sid", "qL","pL","q1","q2","p1","p2", "evalue", "score", "pident", "qcovhsp"))
BlastpOut = filter(BlastpOut, sid %in% na.omit(WorkTree$tip.label)) # Only keep relevant subject (matches) 
BlastpOut = filter(BlastpOut, qid %in% na.omit(IDDT$RID[wna(IDDT$RCR90)])) # Only keep relevant queries - contigs with RdRPs that aren't already affiliated.
setorderv(setDT(BlastpOut),c("qid","score", "pident", "evalue"), order = c(-1,-1,-1,1)) # Sort the hits by alignment stats. 
BlastpOut = unique(BlastpOut, incomparables = FALSE, fromLast = FALSE, by="qid") # Cull the matches to the best hit per query RdRP.
BlastpOut =  Rename1Col(BlastpOut,"qid","RID")

BlastpOut = filter(BlastpOut,pident >= 90) # Remove matches with identity below 90 percent.
BlastpOut = filter(BlastpOut,qcovhsp >= 75) # Remove matches with query coverage below 75 percent.
BlastpOut$AfLvl = p0("Lvl B. - BLASTp match ID >= 90% and qcov >= 75% to megatree leaf:  ", BlastpOut$sid )
BlastpOut =  Rename1Col(BlastpOut,"sid","RCR90")
BlastpOut$node = unlist(NLabel2NID(WorkTree,BlastpOut$RCR90))
BlastpOut = merge(BlastpOut,`dropcols<-`(IDDT,c("RCR90","node","AfLvl")),by="RID",all.x=T,all.y=F)

tmpdf <- filter(IDDT[wna(IDDT$node),], ND %ni% BlastpOut$ND)

ShCl <- SharedCols(IDDT,BlastpOut)
IDDT <- distinct(rbindlist(list(IDDT[wana(IDDT$RCR90),],BlastpOut[,..ShCl],tmpdf)))

#### Level C. - Same RvANI90 cluster as level {A,B} ####
# Create a table with RvANI90 clusters that are associcated with only a single tree leaf.
ANI90tbl = IDDT %>% group_by(RvANI90) %>% summarise(
  RvANI90 = unique(RvANI90),
  RCR90_list = list(unique(RCR90)),
  nodes = list(as.num(na.omit(unique(node)))))
ANI90tbl = ANI90tbl[wana(ANI90tbl$RvANI90),]
ANI90tbl$Nnodes = lengths((ANI90tbl$nodes))
ANI90tbl$Nnodes[wna(ANI90tbl$nodes)] = 0
ANI90tbl <- ANI90tbl[wh(ANI90tbl$Nnodes ==1),]
ANI90tbl$RCR90 <- na.omit(unlist(ANI90tbl$RCR90_list))

W1 <- intersect(wh(IDDT$RvANI90 %in% ANI90tbl$RvANI90),wna(IDDT$RCR90))
tmpdf <- IDDT[W1,]
IDDT <- IDDT[-W1,]

tmpdf = merge(`dropcols<-`(tmpdf,c("RCR90","node","AfLvl")),ANI90tbl[,c("RvANI90","RCR90")] ,by="RvANI90",all.x=T,all.y=F)
tmpdf$AfLvl = p0("Lvl C. - Same RvANI90 cluster as level {A,B}")
tmpdf$node = unlist(NLabel2NID(WorkTree,tmpdf$RCR90))
IDDT <- distinct(rbindlist(list(IDDT,tmpdf),use.names=TRUE))
wna(IDDT$RCR90)

# BlastpOut = merge(BlastpOut,`dropcols<-`(IDDT,c("RCR90","node","AfLvl")),by="RID",all.x=T,all.y=F)
# 
# tmpdf <- filter(IDDT[wna(IDDT$node),], ND %ni% BlastpOut$ND)
# 
# ShCl <- SharedCols(IDDT,BlastpOut)
# IDDT <- distinct(rbindlist(list(IDDT[wana(IDDT$RCR90),],BlastpOut[,..ShCl],tmpdf)))

#### Level D. - best MEGABLAST hit from contig to {A,B,C}; pident >=90% AND (qcov >=75% OR nident >=900 nt) ####
# For this step, we will use several blastn results; the key reason is to avoid limiting the number of matches, or allowing matches to contigs that were eventually filtered off in quality control.
DCB <- fread("../Wolf/Tmp_vs_IDDT3_results.tsv", col.names = c("qid", "sid", "qcov", "nid", "pos", "pid", "qL", "pL",  "q1", "q2", "p1", "p2", "length", "evalue", "score", "gaps"))
AvA <- fread("../Wolf/IDFT_blastn_all_vs_all_results.tsv", col.names =c("qid", "sid", "qcov", "nid", "pos", "pid", "qL", "pL",  "q1", "q2", "p1", "p2", "length", "evalue", "score", "gaps"))
DCB <- distinct(plyr::rbind.fill(DCB,AvA))
cont.unaff <- fread("../Wolf/cont.unaff.0.tab", col.names = c("qid", "sid", "qL", "pL", "q1", "q2", "p1", "p2", "evalue", "score", "pid", "nid", "qcov"))
ShCl <- SharedCols(DCB,cont.unaff)
cont.unaff <- merge(cont.unaff, DCB, by = ShCl, all.x=T,all.y=F)
DCB = distinct(plyr::rbind.fill(DCB[],cont.unaff[]))
DCB <- Rename1Col(DCB, "qid", "ND")
DCB <- filter(DCB, sid %in% IDDT$ND[wana(IDDT$RCR90)]) # We only want matches TO contigs already affiliated.
DCB <- filter(DCB, ND %in%  IDDT$ND[wna(IDDT$RCR90)]) # We only want matches FROM contigs not-already affiliated.

# Filtering to avoid false positive links.
DCB <- DCB[unique(union(wh(DCB$qcov >= 75),wh(DCB$nid >= 900))),]
DCB <- filter(DCB, pid >= 90 ) # Minimum identity of 90 percent.
DCB <- filter(DCB, evalue <= 0.001) # Maximal E-value of 0.001

setorderv(setDT(DCB),c("ND", "score", "nid", "qcov", "evalue"), order = c(-1,-1,-1,-1,1)) # Sort the hits by alignment stats. 
DCB <- unique(DCB, incomparables = FALSE, fromLast = FALSE, by="ND")  # Cull the matches to the best hit per query contig.
DCB = Rename1Col(DCB,"ND","qid")
DCB = Rename1Col(DCB,"sid","ND")

DCB$AfLvl = p0("Lvl D. - best dc-MEGABLAST match with pident >=90% AND (qcov >=75% OR nident >=900nt) to contigs from {A,B,C}:  ", DCB$ND) 
DCB = merge(DCB,IDDT[,c("ND","RCR90")], by="ND", all.x=T, all.y=F)
DCB = DCB[,c("qid","RCR90","AfLvl")]
DCB = Rename1Col(DCB,"qid","ND")
DCB = merge(DCB,`dropcols<-`(IDDT,c("AfLvl","RCR90")),by="ND", all.x=T, all.y=F)


IDDT <- distinct(rbindlist(list(filter(IDDT,ND %ni% DCB$ND),DCB),use.names=TRUE))


fwrite(IDDT,"Novo.IDDT.06022022.tsv",)










