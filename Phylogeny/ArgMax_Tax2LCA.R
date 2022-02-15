# Information  --------------------------------------------------------------------
## Script name: ArgMax_Tax2LCA.R
## Email: uri.neri@gmail.com
##

# Note! 
# This script was includes code that was not necessarily ran in a chronological order.
# Some of the objects were created by other Rscripts in the Phylogentic analysis (such as leaf2tax and cont2leaf).
# As such, it is not advised to attempt running this code 'as is', rather, it is intended to be a disclosure rather than a re-runnable "production" code.
# That said, we would be happy to try and help you should you want to source, modify, reuse, rerun or fork in any way the code.
# (Although, you don't have to tell us etc, this is open source so go ahead and use what you can as much as you want!)
# HTH
# Neri

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
RankDF <- fread("Misc/RankDF.tsv", header = T, sep = "\t")
Rnkcls <- BBmisc::capitalizeStrings(tolower(RankDF$Rank))
RankDF$Rank <- BBmisc::capitalizeStrings(tolower(RankDF$Rank))
#### Read Info ####
IDFT <- Rename1Col(fread0("../Vfin/Metadata/IDFT.14092021.tsv"), "Contig", "ND")
IDFT.fasta <- readDNAStringSet("/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ALL_nuc_3007.fasta")
AllDF <- readRDS(file = "../Annotations/AllDF.07112021.RDS")
WorkTree <- as.phylo(read.tree("./210822/ali.210822m/210822m.77510.09.lab_mod.tre"))

# Misc. ---------------------------------------------------------------
#### 10-17.10.2021 ####

IDFT$Host[which(IDFT$RID == "rv20_v2_2738")] = "Eukaryota; Viridiplantae; Chlorophyta; Ulvophyceae; TCBD clade; Bryopsidales; Bryopsidineae; Bryopsidaceae; Bryopsis; Chloroplast"
IDFT$Host[which(IDFT$RID %in% c("rv20_v2_2876","rv20_v2_4192"))] = "Eukaryota; Viridiplantae; Chlorophyta; Ulvophyceae; TCBD clade; Bryopsidales; Bryopsidineae; Bryopsidaceae; Bryopsis; Mitochondria"

IDFT$Note = gsub(pattern = ", putative bacterial lysis protein encoding",replacement = "",x = IDFT$Note,fixed = T)
IDFT$Note = gsub(pattern = "putative bacterial lysis protein encoding",replacement = "",x = IDFT$Note,fixed = T)




# IDFT = IDFT[-which(IDFT$RID == "Rv4_270391"),]
# IDFT = IDFT[-which(IDFT$ND == "ND_023489"),]

##### Attempt argmaximal LCA of non LCA/LSA'd taxa  #####
TaxCnts = 1:nrow(TaxDT)
TaxNyms = TaxDT$Taxa[TaxCnts]

PramDF = data.frame(Nodetbl[, c("node", "node.label")])
PramDF[, as.char(TaxCnts)] = 0

PramDF[, c("N_c", "M_c", "Kids", "tips", "IsTip")] = NA
PramDF$Kids = Descendants(WorkTree, PramDF$node, type = "all")
PramDF$tips = Descendants(WorkTree, PramDF$node, type = "tips")
PramDF$N_c = lengths(PramDF$tips)
PramDF$IsTip = isTip(WorkTree, PramDF$node)
PramDF$M_c = 0
PramDF$M_c[which(PramDF$IsTip)] = 0.000000000001 # Iffy patch.

Rnkcls = setdiff(intersect(RankDF$Rank, colnames(Nodetbl)), c("Kingdom", "Realm"))
# Rnkcls = intersect(tmpvar$Rank, colnames(Nodetbl))

cntr = 0
pb = txtProgressBar(
  min = 0,
  max = 100,
  initial = 0,
  char = "=",
  width = NA,
  title,
  label,
  style = 3
)

Row2IterOn = which(!(PramDF$IsTip))
Row2IterOn = 1:nrow(PramDF)

# AllTrees = get_subtrees_at_nodes(WorkTree, PramDF$node.label[Row2IterOn])
# AllDists = get_all_distances_to_root(WorkTree)  # Numeric vector of size Ntips+Nnodes, with the i-th element being the distance (cumulative branch length) of the i-th tip or node to the root. Tips are indexed 1,..,Ntips and nodes are indexed (Ntips+1),..,(Ntips+Nnodes).
# {
#   mockParmDF = PramDF
#   sum_
#   PramDF$M_c   get_all_distances_to_root()
#   Nodetbl$DistrFromRoot
# }



# Get nodes M values - Mass index for taxon t (iz) in the tree clade c. (M_t, M_tc and M_c as the sum of branch lengths in the tree pruned for taxon t, the clade c pruned for taxon t, and the whole clade c respectively).
# In PramDF, each taxa represents a column. Additionally, each rank also has a column, which is the sum of all taxa of that rank.
for (i in Row2IterOn) {
  All_Clade_Tips = unlist(PramDF$tips[i])
  w0000 = filter(Nodetbl, node %in% All_Clade_Tips)
  PramDF$M_c[i] = sum(w0000$DistrFromRoot) # Calculates the phylogenetic distance (cumulative branch length) of the root to each tip and node.
  for (ix in Rnkcls) {
    iy = which(TaxNyms %in% (unique(w0000[, ..ix][[1]])))
    if (len(iy) == 0) {
      next
    }  # print("NFF");
    for (iz in iy) {
      w0001 = which(w0000[, ..ix] == TaxNyms[iz])
      PramDF[i, as.char(TaxCnts[iz])] =  sum(w0000$DistrFromRoot[w0001])
      # print(i)
    }
  }
  cntr = cntr + 1
  setTxtProgressBar(pb, as.num(format((
    cntr * 100 / (len(Row2IterOn))
  ), digits = 5)))
}

colnames(PramDF) = c("node",
                     "node.label",
                     TaxNyms,
                     "N_c",
                     "M_c",
                     "Kids",
                     "tips",
                     "IsTip")
PramDF[, Rnkcls] = 0
for (ix in setdiff(Rnkcls, c("Realm", "superkingdom", "Kingdom"))) {
  iy  <-  TaxNyms[which(TaxNyms %in% filter(TaxDT, Rank == ix)$Taxa)]
  PramDF[, ix] <- apply(
    PramDF[, iy],
    MARGIN = 1,
    FUN = function(x)
      sum(x),
    simplify = T
  )
}
saveRDS(PramDF,"PramDF.RDS")
fwrite(sep = '\t',file = "PramDF.tbl",x = PramDF,quote = T,verbose = T)


# Get nodes C values - coverage index for taxon t (iz) in the tree clade c.  C_tc = M_tc / M_t; the coverage index for taxon t in the tree clade c where M_tc is the mass of leaves in clade c belonging to taxon t and M_t is the total mass of leaves belonging to taxon t in the tree.
C_tb.tbl <- `dropcols<-`(PramDF, c(TaxNyms))
C_tb.tbl[, c(TaxNyms)] = 0
for (ix in Rnkcls) {
  iy = TaxNyms[which(TaxNyms %in% filter(TaxDT, Rank == ix)$Taxa)]
  for (iz in iy) {
    C_tb.tbl[, iz] =  PramDF[, iz] /  TaxDT$M_t[which(TaxDT$Taxa == iz)] # PramDF$M_c
  }
}
saveRDS(C_tb.tbl,"C_tb.tbl.RDS")
fwrite(sep = '\t',file = "C_tb.tbl",x = C_tb.tbl,quote = T,verbose = T)
C_tb.tbl[150090, iz] 

# Get nodes P values - "Purity" index; P_tc = M_tc / M_c ; where M_c is the sum of distances to the root of leaves in clade c, and M_tc is the sum of a subset of those leaves, that are annotated as taxa t.
P_tb.tbl <- `dropcols<-`(PramDF, c(TaxNyms))
P_tb.tbl[, c(TaxNyms)] = 0
for (ix in Rnkcls) {
  iy = TaxNyms[which(TaxNyms %in% filter(TaxDT, Rank == ix)$Taxa)]
  for (iz in iy) {
    P_tb.tbl[, iz] = PramDF[, iz] / PramDF[, ix]
    P_tb.tbl[which(is.nan(P_tb.tbl[, iz])), iz] = 0
  }
}
saveRDS(P_tb.tbl,"P_tb.tbl.RDS")
fwrite(sep = '\t',file = "P_tb.tbl",x = P_tb.tbl,quote = T,verbose = T)

PramDF[150090, ix]

# P_tb.tbl[150090, iz]
# 6684
# Get nodes Q values -
Q_tb.tbl <- `dropcols<-`(PramDF, c(TaxNyms))
Q_tb.tbl[, c(TaxNyms)] <- 0
for (ix in Rnkcls) {
  iy = TaxNyms[which(TaxNyms %in% filter(TaxDT, Rank == ix)$Taxa)]
  for (iz in iy) {
    Q_tb.tbl[, iz] = C_tb.tbl[, iz] * P_tb.tbl[, iz]
    Q_tb.tbl[which(is.nan(Q_tb.tbl[, iz])), iz] = 0
  }
}
saveRDS(Q_tb.tbl,"Q_tb.tbl.RDS")
fwrite(sep = '\t',file = "Q_tb.tbl",x = Q_tb.tbl,quote = T,verbose = T)

TaxDT$MaxQtc_Node <- ""
TaxDT[, c("MaxQtc", "MaxQtc_C", "MaxQtc_P")] = 0
Row2IterOn = TaxCnts #TaxDT[, which(!IsMonophyletic_wNF)]

for (i in Row2IterOn) {
  # Shouldn't grow variables like this.
  w01 = which(Q_tb.tbl[, TaxDT$Taxa[i]] == max(Q_tb.tbl[, TaxDT$Taxa[i]]))
  if (len(w01) > 1) {
    print(i)
    w01 = w01[which(C_tb.tbl[w01, TaxDT$Taxa[i]] == max(C_tb.tbl[w01, TaxDT$Taxa[i]]))]
  }
  if (len(w01) > 1) {
    w01 = w01[which(P_tb.tbl[w01, TaxDT$Taxa[i]] == max(P_tb.tbl[w01, TaxDT$Taxa[i]]))]
  }
  if (len(w01) > 1) {
    print(i)
    
    w01 = w01[which.max(Nodetbl$DistrFromRoot[w01])]
  }
  if (length(w01) == 0) {
    next
  }
  TaxDT$MaxQtc_Node[i] = Q_tb.tbl$node.label[w01]
  TaxDT$MaxQtc[i] = Q_tb.tbl[w01, TaxDT$Taxa[i]]
  TaxDT$MaxQtc_P[i] = P_tb.tbl[w01, TaxDT$Taxa[i]]
  TaxDT$MaxQtc_C[i] = C_tb.tbl[w01, TaxDT$Taxa[i]]
}

w1 = which(TaxDT$MaxQtc_Node == "")
len(w1)

View(TaxDT[, c(
  "Taxa",
  "Rank",
  "MaxQtc_Node",
  "MaxQtc",
  "M_t",
  "LCA",
  "LCA_label",
  "MaxQtc_C",
  "MaxQtc_P",
  "IsMonophyletic_wNF"
)])

TaxDT4W = `dropcols<-`(TaxDT, c("UnwantedKids" , "WantedKids"))

TaxDT4W = Rename1Col(dfdt = TaxDT4W,
                     colnm = "LCA",
                     newcolnm = "ABS_LCA")
TaxDT4W = Rename1Col(dfdt = TaxDT4W,
                     colnm = "LCA_label",
                     newcolnm = "ABS_LCA_label")
TaxDT4W = Rename1Col(dfdt = TaxDT4W,
                     colnm = "MaxQtc_Node",
                     newcolnm = "MaxQtc_node.label")
TaxDT4W$MaxQtc_Node = treeio::nodeid(WorkTree, TaxDT4W$MaxQtc_node.label)
TaxDT4W = TaxDT4W[which(TaxDT4W$ABS_LCA_label != "ASD"),]

# View(TaxDT4W)

WriteWolfTbl(TaxDT4W, "TaxDT.19102021.tsv")

write.tree(WorkTree, "Mod.tmp.19102021.tre")
WriteWolfTbl(Nodetbl, "Nodetbl.19102021.tsv")

# Populate Nodetbl with the argmaximal LCAs
NodetblArgMax = data.table(Nodetbl)
for (i in 1:nrow(TaxDT4W)) {
  if(TaxDT4W$MaxQtc[i] <= 0.5){
    print(i)
    next
  }
  tmprank = which(colnames(NodetblArgMax) == TaxDT4W[i, Rank])
  tmpkids = unique(unlist(GetKids(WorkTree, TaxDT4W$MaxQtc_Node[i], type = "all")))
  tmplabs = NID2NLabel(WorkTree, tmpkids)
  NodetblArgMax[node.label %in% tmplabs, tmprank] = as.char(TaxDT4W[i, Taxa])
}

tmpcold = c("node.label", Rnkcls)

NodetblComp = merge(
  Nodetbl,
  NodetblArgMax[, ..tmpcold],
  by = "node.label",
  all.x = T,
  all.y = F
)
WriteWolfTbl(NodetblComp, "NodetblComp.19102021.tsv")
WriteXlsx(NodetblComp, "NodetblComp.19102021.xlsx")

whd(NodetblArgMax$node.label, returnValue = T)

fwrite(sep = '\t',file = "NodetblArgMax.tbl",x = NodetblArgMax,quote = T,verbose = T)
fwrite(sep = '\t',file = "IDDT6.tbl",x = IDDT6,quote = T,verbose = T)

{
  sh1 = SharedCols(IDDT6,IDFT)
  sh2 = SharedCols(IDDT6,NodetblArgMax)
  ttt1 = setdiff(c(sh2,sh1),RankDF$Rank)
  IDDT7 = merge(IDDT6[,ttt1],IDFT[,c("Host","ND","Note","Source","Genetic_Code","Spid","Type.hit","Hit.s.","Full_name","Full_name2","RBS")],by = "ND",all.x=T, all.y=F)
  IDDT7$Host[which(IDDT7$NC90 %in% filter(NodetblArgMax,Class ==  "Leviviricetes")$NC90)] = "Bacteria"
  IDDT7$Host[which(IDDT7$NC90 %in% filter(NodetblArgMax,Class ==  "Vidaverviricetes")$NC90)] = "Bacteria"
  w1 = which(IDDT7$ND %in% Putative_Lysis_Encoding)
  IDDT7$Note[w1] = p0("Putative_Lysis_Encoding, ", IDDT7$Note[w1])
  fwrite(sep = '\t',file = "IDDT7.tbl",x = IDDT7,quote = T,verbose = T)
  
}
####### (Deprecated!) Revert node affiliation for immutable taxa/known clades. (Deprecated!) ####### 
# NodetblArgMax_NoRev = NodetblArgMax
# setDF(Nodetbl)
# setDF(NodetblArgMax)
# AllRnkcls = Rnkcls
# for (ix in AllRnkcls) {
#   w01 = which(Nodetbl[, ix] != "NF")
#   if (length(w01) != 0) {
#     NodetblArgMax[w01, ix] = Nodetbl[w01, ix]
#   }
# }
# WriteWolfTbl(`dropcols<-`(NodetblArgMax, c("sseqid", "taxid")), "NodetblArgMax.210822m.tsv")
# WriteWolfTbl(`dropcols<-`(NodetblArgMax_NoRev, c("sseqid", "taxid")),"NodetblArgMax_NoRev.210822m.tsv")
# write.tree(WorkTree, "Mod.tmp.210822m.tre")
# WriteWolfTbl(`dropcols<-`(Nodetbl, c("sseqid", "taxid")), "Nodetbl.210822m.tsv")
