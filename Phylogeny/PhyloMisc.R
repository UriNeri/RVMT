# Information  --------------------------------------------------------------------
## Script name: PhyloMisc.R
## Author: Uri Neri
## Email: uri.neri@gmail.com
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


# Identify "misbehaving" taxa  --------------------------------------------------------------------
# TmpRnkChar = SharedCols(VMR_MSL36,NodetblArgMax)
TmpRnkChar = c("Phylum","Class","Order","Family","Genus")
LineageTaxnmy = distinct(VMR_MSL36[,TmpRnkChar])
TmpNtwAm = data.frame(NodetblArgMax)
Check_TaxDT4W = filter(TaxDT4W,MaxQtc > 0.5)
Check_TaxDT4W = filter(Check_TaxDT4W,Rank %in% colnames(LineageTaxnmy))

Check_TaxDT4W$IsModLin = F
CDTF = Check_TaxDT4W[,c("Taxa","Rank","IsModLin")]
CDTF[,c("Old_Lineage","New_lineage")] = 0
for (i in 1:nrow(Check_TaxDT4W)) {
  tmprank = which(colnames(TmpNtwAm) == Check_TaxDT4W[i, Rank])
  ParentRank = TmpRnkChar[which(TmpRnkChar == Check_TaxDT4W[i, Rank])-1]
  tmpParentRank = which(colnames(TmpNtwAm) == ParentRank)
  
  w1 = which(TmpNtwAm[,tmprank] == Check_TaxDT4W[i, Taxa])
  ObservedParent = unique(TmpNtwAm[w1,tmpParentRank])
  
  w2 = which(LineageTaxnmy[,Check_TaxDT4W[i, Rank]] == Check_TaxDT4W[i, Taxa])
  ExpectedParent = unique(LineageTaxnmy[w2,ParentRank] )
  
  w4 = wana(unique(ObservedParent))
  w5 = wana(unique(ExpectedParent))
  
  if((len(w4) == 0) & (len(w5 != 0 ))){
    print(p0("For ", Check_TaxDT4W[i, Taxa], "  observed parent rank is NA, but expected is ", toString(unique(ExpectedParent))))
    Check_TaxDT4W$IsModLin[i] = T
    
    # Old_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(LineageTaxnmy[w2,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    # New_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(TmpNtwAm[w1,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    Old_Lineage = toString(p0(ExpectedParent,";",Check_TaxDT4W[i, Taxa]))
    New_Lineage = toString(p0(ObservedParent,";",Check_TaxDT4W[i, Taxa]))
    
    CDTF$Old_Lineage[i] = Old_Lineage
    CDTF$New_Lineage[i] = New_Lineage
    next
  }
  if((len(w5) == 0) & (len(w4 != 0 ))){
    print(p0("For ", Check_TaxDT4W[i, Taxa], "  expected parent rank is NA, but observed is ", toString(unique(ObservedParent))))
    Check_TaxDT4W$IsModLin[i] = T
    # Old_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(LineageTaxnmy[w2,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    # New_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(TmpNtwAm[w1,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    Old_Lineage = toString(p0(ExpectedParent,";",Check_TaxDT4W[i, Taxa]))
    New_Lineage = toString(p0(ObservedParent,";",Check_TaxDT4W[i, Taxa]))
    
    CDTF$Old_Lineage[i] = Old_Lineage
    CDTF$New_Lineage[i] = New_Lineage
    next
  }
  if((len(w5) == 0) & (len(w4) == 0 )){
    # print(p0("For ", Check_TaxDT4W[i, Taxa], "  both expected and observed parent ranks are NA"))
    Check_TaxDT4W$IsModLin[i] = F
    # Old_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(LineageTaxnmy[w2,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    # New_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(TmpNtwAm[w1,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    Old_Lineage = toString(p0(ExpectedParent,";",Check_TaxDT4W[i, Taxa]))
    New_Lineage = toString(p0(ObservedParent,";",Check_TaxDT4W[i, Taxa]))
    
    CDTF$Old_Lineage[i] = Old_Lineage
    CDTF$New_Lineage[i] = New_Lineage
    next
  }
  
  if (len(w4) > 1 ){
    print(p0("For ", Check_TaxDT4W[i, Taxa], "  expected to be under ", ExpectedParent, " but is under multiple observed taxa - ", toString(unique(ObservedParent))))
    Check_TaxDT4W$IsModLin[i] = T
    # Old_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(LineageTaxnmy[w2,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    # New_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(TmpNtwAm[w1,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    Old_Lineage = toString(p0(ExpectedParent,";",Check_TaxDT4W[i, Taxa]))
    New_Lineage = toString(p0(ObservedParent,";",Check_TaxDT4W[i, Taxa]))
    
    CDTF$Old_Lineage[i] = Old_Lineage
    CDTF$New_Lineage[i] = New_Lineage
    next
  }
  
  if (unique(na.omit(ObservedParent)) != unique(na.omit(ExpectedParent))){
    print(p0("For ", Check_TaxDT4W[i, Taxa], "  expected to be under ", ExpectedParent, " but is under ", toString(unique(ObservedParent))))
    Check_TaxDT4W$IsModLin[i] = T
    # Old_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(LineageTaxnmy[w2,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    # New_Lineage = toString(BBmisc::rowSapply(df = data.frame(distinct(TmpNtwAm[w1,TmpRnkChar])),fun = function(x) toString(c(x,"|||"))))
    Old_Lineage = toString(p0(ExpectedParent,";",Check_TaxDT4W[i, Taxa]))
    New_Lineage = toString(p0(ObservedParent,";",Check_TaxDT4W[i, Taxa]))
    
    CDTF$Old_Lineage[i] = Old_Lineage
    CDTF$New_Lineage[i] = New_Lineage
    next
  }
} 
CDTF$IsModLin = Check_TaxDT4W$IsModLin
CDTF = RemoveEmptyStrRow(CDTF)
CDTF = CDTF[CDTF$IsModLin,]
ObservedDF <- plyr::count(distinct(TmpNtwAm[,c("Phylum","Class","Order","Family","Genus")]) )
ExpectedDF <- plyr::count(distinct(LineageTaxnmy[,c("Phylum","Class","Order","Family","Genus")]) )
saveRDS(CDTF,"CDTF.RDS")
fwrite(sep = '\t',file = "CDTF.tbl",x = CDTF,quote = T,verbose = T)
WriteWolfTbl(CDTF,"CDTF.tsv")
WriteWolfTbl(ObservedDF,"ObservedDF.tsv")
WriteWolfTbl(ExpectedDF,"ExpectedDF.tsv")


# Inspect to see if any sub-group label occurred in more than one parent group:
tmpdf = plyr::count(distinct(TmpNtwAm[,c("Class","Phylum")]),vars ="Class" )
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Class","Order")]),vars ="Order" )
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Family","Order")]),vars ="Family" )
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Family","Genus")]),vars ="Genus" )


# Autoname unaffiliated clades  --------------------------------------------------------------------
NodetblArgMax.wAN = NodetblArgMax # Nodetbl argmaxd with AutoNamed
TreeTrav = get_tree_traversal_root_to_tips(WorkTree, include_tips = T)$queue
NrootTrav = TreeTrav[2:length(TreeTrav)]
# NrootTrav = TreeTrav
NewTaxDT = data.frame(
  "Taxa" = rep("", 10000),
  "C_tc" = rep(0, 10000),
  "LCA" = rep(0, 10000),
  "LCA_label" = rep("", 10000),
  "Rank" = rep("", 10000),
  "Cons" = rep(0, 10000)
) # , UnwatedKids = vector('list',10000))
NewTaxDT$WantedKids = vector('list', nrow(NewTaxDT))   # List place holder 1
NewTaxDT$UnwantedKids = vector('list', nrow(NewTaxDT)) # List place holder 2
NewTaxDT$Kids = vector('list', nrow(NewTaxDT)) # List place holder 3


RTlca = get_mrca_of_set(WorkTree, filter(taxiddf_lin1, Kingdom %in% "RTs")$node.label)
RTsNodes = c(GetKids(WorkTree, RTlca, type = "all"), RTlca)
TmC = 1
# NodesDone = c(77568,77569)
NodesDone = c(-1)
UpperRankNodesD = 0
# for (ix in Rnkcols[-c(1,6)]){
for (ix in c("Phylum", "Subphylum", "Class", "Order", "Family")) {
  for (iy in NrootTrav) {
    if (iy %in% NodesDone) {
      next
    }
    if (iy %in% RTsNodes) {
      next
    }
    if (iy %in% UpperRankNodesD) {
      next
    }
    Nlb = nodelab(WorkTree, iy) # Node label
    Nyd = which(NodetblArgMax.wAN$node.label == Nlb) # Node index in the Nodetbl
    if (NodetblArgMax.wAN[Nyd, ix] != "NF") {
      next
    }
    NewTaxDT$LCA[TmC] = iy
    NewTaxDT$LCA_label[TmC] = Nlb
    NewTaxDT$Rank[TmC] = ix
    NewTaxDT$Kids[TmC] = list(NID2NLabel(WorkTree, GetKids(WorkTree, iy, type =
                                                             "all")))
    w00 = which(NodetblArgMax.wAN$node.label %in%  unlist(NewTaxDT$Kids[TmC]))
    w01 = NodetblArgMax.wAN[w00,]
    NewTaxDT$UnwantedKids[TmC] = list(w01$node.label[which(w01[, ix] != "NF")])
    NewTaxDT$WantedKids[TmC] = list(w01$node.label[which(w01[, ix] == "NF")])
    NewTaxDT$Cons[TmC] = lengths(NewTaxDT$WantedKids[TmC]) / lengths(NewTaxDT$Kids[TmC])
    Temperere = w01$DistrFromRoot[which(w01[, ix] != "NF")]
    Temperere2 = w01$DistrFromRoot[which(w01[, ix] == "NF")]
    if (len(Temperere) == 0) {
      NewTaxDT$C_tc[TmC] = 1
    }
    if (len(Temperere) != 0) {
      NewTaxDT$C_tc[TmC] = sum(Temperere2)  / (sum(Temperere) + sum(Temperere2))
    }
    w012 = which(NodetblArgMax.wAN$node.label %in% unlist(NewTaxDT$WantedKids[TmC]))
    NewTaxDT$Taxa[TmC] = p0("C.", iy, ".R.", ix)
    NodesDone = unique(c(NodesDone, iy))
    if (NewTaxDT$C_tc[TmC] < 0.5) {
      NodesDone = unique(c(NodesDone, NLabel2NID(
        WorkTree, unlist(NewTaxDT$UnwantedKids[TmC])
      )))
      NewTaxDT[TmC,] = NewTaxDT[TmC + 1,]
      next
    }
    NodesDone = unique(c(NodesDone, c(iy, GetKids(
      WorkTree, iy, type = "all"
    ))))
    NodetblArgMax.wAN[c(Nyd, w012), ix] = NewTaxDT$Taxa[TmC]
    TmC = TmC + 1
    # NrootTrav = setdiff(NrootTrav,GetKids(WorkTree,iy,type="all"))
  }
  NodesDone = c(-1)
  UpperRankNodesD = unique(c(UpperRankNodesD, NewTaxDT$LCA))
  # NrootTrav = TreeTrav[2:length(TreeTrav)]
  print(p0(iy, "      -       ", TmC))
}
NewTaxDT = distinct(NewTaxDT)
NewTaxDT4W = `dropcols<-`(NewTaxDT, c("UnwantedKids" , "WantedKids", "Kids"))
NewTaxDT4W = Rename1Col(dfdt = NewTaxDT4W,
                        colnm = "LCA",
                        newcolnm = "node")
NewTaxDT4W = Rename1Col(dfdt = NewTaxDT4W,
                        colnm = "LCA_label",
                        newcolnm = "label")
NewTaxDT4W$DistFromRoot = get_all_distances_to_root(WorkTree, as_edge_count = F)[NLabel2NID(WorkTree, NewTaxDT4W$label)]
WriteWolfTbl(NewTaxDT4W, "NewTaxDT.210822m.tsv")
NodetblArgMax.wAN$sseqid = NULL
NodetblArgMax.wAN$taxid = NULL

WriteWolfTbl(NodetblArgMax.wAN, "NodetblArgMax.wAN.210822m.tsv")
NewTaxDT = NewTaxDT[which(NewTaxDT$Cons > 0),]


View(NewTaxDT[, c("Taxa", "C_tc", "LCA", "LCA_label",  "Rank", "Cons")])
gc()

# #Plots  --------------------------------------------------------------------
library(ggtree)
library(ggrepel)
# library(plotly)
# library(ggplotify)
library(ggtreeExtra)
library(ggimage)
# library(treeio)
# library(tidytree)
library(ggstar)
# library(ggplot2)
library(ggnewscale)

###### Make ggtrees ###### 
# ggWorkTreeCirc = ggtree(WorkTree, layout = "circular", yscale = "none", ladderize = T)
# saveRDS(ggWorkTreeCirc, "ggWorkTreeCirc.phylo.tre.RDS")
ggWorkTreeCirc = readRDS("ggWorkTreeCirc.phylo.tre.RDS")

# ggWorkTreeApe = ggtree(WorkTree, layout = "ape", yscale = "none", ladderize = T)
# saveRDS(ggWorkTreeApe, "ggWorkTreeApe.phylo.tre.RDS")
ggWorkTreeApe = readRDS("ggWorkTreeApe.phylo.tre.RDS")

# ggWorkTreeRect = ggtree(WorkTree, layout = "rectangular", yscale = "none", ladderize = T)
# saveRDS(ggWorkTreeRect, "ggWorkTreeRect.phylo.tre.RDS")
ggWorkTreeRect = readRDS("ggWorkTreeRect.phylo.tre.RDS")


ggWorkTreeCircUltra <- readRDS("ggWorkTreeCircUltra.RDS.2022-01-12.17-25-13.RDS")
###### Prepare node annotation data ###### 
IDFT$BinID <- p0(
  str_split_fixed(string = IDFT$ND, pattern = fixed("_"),n = 3)[,1],
  "_",
  str_split_fixed(string = IDFT$ND, pattern = fixed("_"),n = 3)[,2])
IDFT[,MlenBin := sum(.SD),.SDcols="Length",by = (BinID)]


tmpAllDFs <- distinct(data.frame(rbindlist(list(AllDF,AllDF6Frx),fill=T))[,SharedCols(AllDF,AllDF6Frx)])
# Movement protein info:
w1 = grep(pattern = "MP_",x = tmpAllDFs$New_Name,fixed=T)
Putative_Movement_Encoding = unique(tmpAllDFs$ND[w1])
Putative_Movement_Encoding_Nodes = unique(na.omit(IDFT$node[which(IDFT$ND %in% Putative_Movement_Encoding)]))
w1 = intersect(wna(IDFT$Host),which(IDFT$node %in% Putative_Movement_Encoding_Nodes))
IDFT$Host[w1] <- "Eukaryota"
IDFT$Host_Evidence[w1] <- "Host related domain - Movement protein"


# Lysis info:
w1 = grep(pattern = "Lysis_sgl|Lysin_|LYZF1",x = tmpAllDFs$New_Name)
w2 = filter(tmpAllDFs[w1,], evalue < 0.005)
Putative_Lysis_Encoding = unique(w2$ND)
Putative_Lysis_Encoding_Nodes = unique(na.omit(IDFT$node[which(IDFT$ND %in% Putative_Lysis_Encoding)]))
Putative_Lysis_Encoding_RID = unique(na.omit(IDFT$node.label[which(IDFT$ND %in% Putative_Lysis_Encoding)]))
w1 = intersect(wna(IDFT$Host),which(IDFT$node %in% Putative_Lysis_Encoding_Nodes))
IDFT$Host[w1] <- "Bacteria"
IDFT$Host_Evidence[w1] <- "Putative Host related domain - bacteriolytic protein"

# CRISPR info:
w1 = union(wana(IDFT$Hit.s.),wana(IDFT$Type.hit))
Putative_CRISPR_Matching_Nodes = unique(na.omit(IDFT$node[w1]))
IDFT$Host_Evidence[wh(IDFT$Host_Evidence == "Putative CRISPR spacer match")] = NA
IDFT$Host[wh(IDFT$Host_Evidence == "Putative CRISPR spacer match")] = NA
w1 = intersect(wna(IDFT$Host_Evidence),wh(IDFT$node %in% Putative_CRISPR_Matching_Nodes))
IDFT$Host[w1] <- "Bacteria"
IDFT$Host_Evidence[w1] <- "Putative CRISPR spacer match"

tmpIDFT <- IDFT
tmpIDFT$Novel <- T
tmpIDFT$Novel[grep(pattern = "al2_",x = tmpIDFT$RID,fixed = T)] = F
tmpIDFT$Novel[grep(pattern = "ap2_",x = tmpIDFT$RID,fixed = T)] = F
tmpIDFT$Novel[grep(pattern = "rv20_",x = tmpIDFT$RID,fixed = T)] <-   F
tmpIDFT$Novel[grep(pattern = "ya2_",x = tmpIDFT$RID,fixed = T)] = F


info <- distinct(tmpIDFT[wana(tmpIDFT$node),] %>% group_by(node) %>% summarise(
  node = node,
  Mlen = max(MlenBin),
  CRHit = toString(unique(na.omit(Type.hit))),
  Host = toString(unique(Host)),
  RBS = toString((RBS)),
  Phylum = unique(Phylum),
  Class = unique(Class),
  Order = unique(Order),
  Family = unique(Family),
  Genus = toString(unique(Genus)),
  # Novel = ((len(wh(Novel))) / len(Novel)  ),
  Novel = AllTrue(Novel),
  Genetic_Codes = toString(unique(Genetic_Code))))

info$node.label <- unlist(NID2NLabel(WorkTree,TmpNodeTblxIDDT$node))

info$Note <- NA
info$Note[info$node %in% Putative_Lysis_Encoding_Nodes]  <-  LysisSign
info$Note[info$node %in% Putative_CRISPR_Matching_Nodes] <-  CRISPRSign

saveRDS(info,"../Wolf/ggtreeinfo.12012022.RDS")

setDF(info)

infoInternals <- data.frame("Taxa"="","rank"="","node"="","node.label"="","maxlen"=0)[0,]
for (rnkcl in c("Phylum","Class","Order","Family","Genus")) {
  AllRnkTaxs <- data.frame(Taxa=unname(unique(info[,rnkcl])),rank=rnkcl)
  infoInternals <- rbind.fill(infoInternals,AllRnkTaxs)
}

setDT(infoInternals)

for (ix in 1:nrow(infoInternals)) {
  tips_ix <- unique(info$node[wh(info[,as.char(infoInternals$rank[ix])] == infoInternals$Taxa[ix])])
  tips_iy <- unlist(NID2NLabel(WorkTree,tips_ix))
  Mlen <- max(IDFT$MlenBin[wh(IDFT$node.label %in% tips_iy)])
  Taxlca <- get_mrca_of_set(WorkTree,tips_ix)
  infoInternals[ix,node:=Taxlca]
  infoInternals[ix,maxlen:=  Mlen]
}
infoInternals[,node:=  as.num(node)]
infoInternals[,node.label:=  unlist(NID2NLabel(WorkTree,node))]

infoInternals2 <- infoInternals[wh(infoInternals$rank != "Genus"),]
infoInternals3 <- infoInternals[grep(pattern = "genPartiti",x = infoInternals$Taxa),]
infoInternals <- rbind(infoInternals2,infoInternals3)
rm(infoInternals2,infoInternals3)

RTsLCA <- get_mrca_of_set(WorkTree,grep(x = WorkTree$tip.label,pattern = "rt",value = T))
infoInternals <- rbind(infoInternals,data.frame("node" = RTsLCA, "node.label" =NID2NLabel(WorkTree,RTsLCA),"rank"="Phylum","Taxa"="Reverse_Transcriptases","maxlen"=0))
w1 <- intersect(wh(infoInternals$rank != "Class"),grep(pattern = fixed(".base-"),x =infoInternals$Taxa,fixed = T))
infoInternals$NymN <- infoInternals$Taxa
infoInternals$NymN[w1] <- str_split_fixed(string = infoInternals$Taxa[w1],pattern = fixed(".base-"),n=2)[,1]
ToKeep <- intersect(wh(info$Note == "X"),unique(c(wh(info$Phylum %in% c("Lenarviricota","p.0002.base-Kitrino")),wh(info$Class == "Duplopiviricetes"))))

Nodetbl <- Tree2Edgetbl(WorkTree)
info2 <- rbind.fill(data.frame(info),data.frame(filter(Nodetbl,!(node.label %in% info$node.label)))[,c("node.label","node")])
infoInternals <- filter(infoInternals,!(Taxa  %in% c("Alvernaviridae", "Barnaviridae", "Carmotetraviridae", "Quadriviridae", "Sunviridae")))

setDT(info2)
for(ix in 1:nrow(infoInternals)){
  Allklids <- unlist(GetKids(WorkTree,node = as.num(infoInternals$node[ix]),type="all",Return_Labels = T))
  iy = wh(info2$node.label %in% Allklids)
  info2[iy, as.char(infoInternals$rank[ix]) :=  infoInternals$Taxa[ix]]
}
setDF(info2)
for(iy in 1:ncol(info2)){
  w1 <- wh(unlist(unname(info2[,iy])) == "NA")
  if(len(w1) !=0){
    info2[w1,iy] = NA
  }
}

infoInternals$DistFromRoot = get_all_distances_to_root(WorkTree, as_edge_count = F)[NLabel2NID(WorkTree, infoInternals$node.label)]
max(infoInternals$DistFromRoot)

info2$DistFromRoot = get_all_distances_to_root(WorkTree, as_edge_count = F)[NLabel2NID(WorkTree, info2$node.label)]
wna(info2$Novel)

AltSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/Alt_Sign.png")
PermuteSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/Permute_Sign.png")
CRISPRSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/BlueStar.png")
LysisSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/RedStar.png")

AltSets <- fread("Alt_Clades.tab",sep = "\t")
AltSets$leaflist = vector('list',c(nrow(AltSets)))   # List place holder 1
AltSets$leaflist <- apply(AltSets,MARGIN = 1,FUN = function(x) as.char(unlist(str_split_fixed(string = x["V2"],pattern = ",",n = Inf))))
AltSets$node <- unlist(apply(AltSets,MARGIN = 1,FUN = function(x) get_mrca_of_set(WorkTree,unlist(x["leaflist"]))))
AltSets$node.label <- NID2NLabel(WorkTree,AltSets$node)

PermutSets <- fread("Perm_Clades.txt",sep = "\t")
PermutSets$leaflist = vector('list',c(nrow(PermutSets)))   # List place holder 1
PermutSets$leaflist <- apply(PermutSets,MARGIN = 1,FUN = function(x) as.char(unlist(str_split_fixed(string = x["V2"],pattern = ",",n = Inf))))
PermutSets$node <- unlist(apply(PermutSets,MARGIN = 1,FUN = function(x) get_mrca_of_set(WorkTree,unlist(x["leaflist"]))))
PermutSets$node.label <- NID2NLabel(WorkTree,PermutSets$node)

info2$shapo <- NA
info2$shapo[wh(info2$node.label %in% AltSets$node.label)] <- "Alternative_Genetic_Code"
info2$shapo[wh(info2$node.label %in% PermutSets$node.label)] <- "Permuted_RdRP"

info2$Novel[wh(info2$Phylum == "Reverse_Transcriptases")] <- F
info2$Mlen[wh(info2$Phylum == "Reverse_Transcriptases")] <- 0


setDT(info2)
info2$AnyKnown <- !(info2$Novel)
info2$AnyKnown[wh(info2$Phylum == "Reverse_Transcriptases")] <- T
info2$NovelProp <- 0
info2$Total_weight <- 0
info2$Novel_weight = 0

info2$AllKids = GetKids(WorkTree,info2$node,type="tips")
info2$Nkids <- lengths(info2$AllKids)
info2$NkidsProp <- info2$Nkids/Ntip(WorkTree)

for(ix in 1:nrow(info2)){
  iy <- wh(info2$node %in% unlist(info2$AllKids[ix]))
  info2[ix, Total_weight:= sum(info2$DistFromRoot[iy])]
  info2[ix, AnyKnown:= AnyFalse(info2$Novel[iy])]
  iy <- wh(info2$Novel[iy])
  if(len(iy)==0){
    next
  }
  info2[ix,Novel_weight:= sum(info2$DistFromRoot[iy])]
  print(ix)
}

info2$NovelProp = info2$Novel_weight / info2$Total_weight
info2 <- distinct(info2)
saveRDS(info2,"../Wolf/ggtreeinfo.12012022.RDS")
# 
# WorkTree2$node.label = WorkTree$node.label
# ggWorkTreeCircUltra <- ggtree(WorkTree2,layout = "circular",  branch.length = "none")


###### MegaTree in Slanted/Fan layout  - Color /+ label by taxonomic rank ###### 
library(ggtree)
library(ragg)

ggWorkTree = readRDS("ggWorkTreeCircUltra.RDS.2022-01-12.17-25-13.RDS")

colorcols <- c("Lenarviricota" = "#0072B2", "Lenar" = "#0072B2",
               "Pisuviricota" = "#009E73","Pisu" = "#009E73",
               "Kitrinoviricota" = "#F0E442","Kitrino" = "#F0E442",
               "Duplornaviricota" = "#8E8C62","Duplorna" = "#8E8C62",
               "Negarnaviricota" = "#CC79A7","Negarna" = "#CC79A7",
               "p.0002.base-Kitrino" = "#DC143C","p.0002" = "#DC143C","<-----------p.0002----|" = "#DC143C",
               "Phylum" = "#DA8D7A","Family" = "#DAD97A", "Genus" = "#7ADA89",  "Order" = "#7ABADA",  "Class"= "#A17ADA",
               "Alt" = "Red", "Standard" = "Blue",
               "Alternative_Genetic_Code" = 17, "Permuted_RdRP"=16)


{
  FakeTree<- ggtree(rtree(30),layout = "circular")
  FakeTree$data$Mlen <- 0
  FakeTree$data$Mlen[FakeTree$data$isTip] <- sample(x=c(1:40),size = 30)
  FakeTree$data$Phylum <- NA
  FakeTree$data$Phylum[sample(wh(FakeTree$data$isTip),30)] <- sample(names(colorcols),size = 30,replace = T)
  # scale_colour_manual(values = colorcols)
  # ito7cols
  # 
  # FakeTree1 <- FakeTree +scale_fill_manual(values = ito7cols)+ #new_scale_fill()+  #scale_colour_manual(values = colorcols) + 
  #   geom_fruit(geom=geom_bar,stat="identity",
  #              mapping=aes(x=Mlen,fill=Phylum),
  #              pwidth=0.5,orientation="y",
  #              grid.params=list(vline =F),
  #              axis.params=list(title = "ito7",
  #                               axis       = "x",
  #                               title.height = 0.001,
  #                               text.size = 2,
  #                               nbreak     = 10))#+scale_colour_manual(values = colorcols)
  # 
  
  FakeTree1 <- FakeTree +scale_fill_manual(values = colorcols)+ #new_scale_fill()+  #scale_colour_manual(values = colorcols) + 
    geom_fruit(geom=geom_bar,stat="identity",
               mapping=aes(x=Mlen,fill=Phylum),
               pwidth=0.5,orientation="y",
               grid.params=list(vline =F),
               axis.params=list(title = "N6Col",
                                axis       = "x",
                                title.height = 0.001,
                                text.size = 2,
                                nbreak     = 10))#+scale_colour_manual(values = colorcols)
  
  
  FakeTree1$data$Code = NA
  FakeTree1$data$Code[wana(FakeTree$data$label)] = sample(c("Standard","Alt"),30,replace = T)
  FakeTree1$data$Code = NULL
  
  FakeTree3 <-  FakeTree1 + geom_fruit(data = data.frame("tlabel"=FakeTree$data$label[wana(FakeTree$data$label)],"Code" = sample(c("Standard","Alt"),30,replace = T)),
                                       geom=geom_point, mapping=aes(y=tlabel,colour=Code),
                                       pwidth=1115,stat="identity",
                                       grid.params=list(vline =F),
                                       axis.params=list(title = "Alt. Genetic Code",
                                                        axis       = "x"))# +
  # theme(axis.title = element_text(hjust = -5,vjust = -4))
  FakeTree3
  FakeTree4 <- FakeTree3
  FakeTree4$data$shapo <- NA
  FakeTree4$data$shapo[37:42]=sample(c("Alt","Perm"),size = 6,replace = T)
  FakeTree5 <- FakeTree4 + geom_nodepoint(aes(size=(6+(NkidsProp*24)),shape=shapo), na.rm = T,color = "purple", fill='purple')
  FakeTree5
  
  
  
  FakeTree3 <- p12 +   geom_nodepoint(mapping = aes(subset = (!is.na(shapo)),shape=shapo,size = (8+(NkidsProp*24)),fill=shapo),)#,color = "purple", fill='purple') #+
  
  
  
  
  # cvd_emulator(FakeTree2)
  p123=ggarrange(plots = list("ito7"=FakeTree1,"Neri6"=FakeTree2,"cb1"=cvd_grid(FakeTree1),"cb2"=cvd_grid(FakeTree2)),nrow = 2)
  QuickTreePDF()
  p123
  dev.off()
  
  
  
}

# p2 <- ggWorkTree
p2 <- ggWorkTree
# p2 <- ggWorkTree #+scale_color_gradient(aes(color = AnyKnown),low="#565656", high="#22ace7")# + #ggtree(tr = WorkTree, layout = "circular", ladderize = T, mapping = aes(color = AnyKnown))#,low="#565656", high="#22ace7")))
p4 <- p2

p4$data <- merge((p4$data),info2,by = "node",all.x = T,all.y = F)
p4$data <- merge((p4$data),`dropcols<-`(infoInternals,c("node.label","DistFromRoot","Lysis","CRISPR","Nkids","NkidsProp")),by = "node",all.x = T,all.y = F)

p4$data$Nkids <- lengths((GetKids(WorkTree,p4$data$node,type="tips")))
p4$data$NkidsProp <- p4$data$Nkids/Ntip(WorkTree)

p4$data <- distinct(merge(p4$data,data.frame(tmpdfp8)[,c("node","Weight_Novel","AnyKnown")],by="node",all.x=T,all.y=F))
p4$data$Note[Putative_Lysis_Matching_Nodes] = LysisSign
p4$data$Note[wh(p4$data$node %in% Putative_CRISPR_Matching_Nodes)] = CRISPRSign

w1 <- wh(p4$data$isTip)

old <- intersect(w1,wh(p4$data$AnyKnown==1))
novel <- setdiff(w1,wh(p4$data$AnyKnown==1))

p4$data$AnyKnown[novel] <- F
p4$data$AnyKnown[old] <- T
p4$data$AnyKnown <- as.num(p4$data$AnyKnown)
p4$data$Weight_Novel[intersect(w1,wh(p4$data$AnyKnown==1))] = 0 
p4$data$Weight_Novel[setdiff(w1,wh(p4$data$AnyKnown==1))] = 1 #p8$data$DistFromRoot.x[setdiff(w1,wh(p8$data$AnyKnown==1))]

p4$data$Note[wh(p4$data$Note == "X")] = LysisSign
p4$data$Note[wh(p4$data$node %in% Putative_CRISPR_Matching_Nodes)] = CRISPRSign

p9 <- p8 + geom_tree(mapping = aes(color = AnyKnown)) + scale_color_gradient(low="#565656", high="#22ace7")# +
# p9 <- p8 + scale_color_gradient(aes(color = AnyKnown),low="#565656", high="#22ace7")# +




p3 = p4 

# 
# p7 <- p3 +scale_fill_manual(values = colorcols)+
#   geom_fruit(geom=geom_bar,mapping=aes(x=Mlen,fill=Phylum),
#              pwidth=0.5,orientation="y",stat="identity",
#              grid.params=list(vline =F),
#              axis.params=list(title = "Maximal observed Genome Length [nt]",
#                               axis       = "x",
#                               title.height = 0.001,
#                               text.size = 3,
#                               hjust      = -0.2,
#                               vjust      = +0.1,
#                               nbreak     = 10,
#                               line.size = 0.4)) +
#   theme(axis.title = element_text(hjust = -5,vjust = -4))
# p7[["layers"]][[7]][["data"]][["y"]] <- 190
## View(p7$data)
p7 <- p3
p8 <- p7

# 
Nodes2Iter <- wna(p8$data$Novel)
novel <- wh(p8$data$Novel)
old <- wh(!(p8$data$Novel))
p8$data$Novel_weight <- 0
p8$data$Novel_weight[novel] <- 1
p8$data$AnyKnown <- T
p8$data$AnyKnown <- !p8$data$Novel
tmpdfp8 <- data.frame(p8$data[Nodes2Iter,c("node","label","Novel_weight","DistFromRoot.x","AnyKnown")])
tmpdfp8$Total_weight <- 0
tmpdfp8$AllKids = GetKids(WorkTree,tmpdfp8$node,type="tips")
setDT(tmpdfp8)
for(ix in 1:nrow(tmpdfp8)){
  iy <- wh(p8$data$node %in% unlist(tmpdfp8$AllKids[ix]))
  tmpdfp8[ix,Total_weight:= sum(p8$data$DistFromRoot.x[iy])]
  tmpdfp8[ix,AnyKnown:= AnyFalse(p8$data$Novel[iy])]
  iy <- wh(p8$data$Novel[iy] ==1)
  if(len(iy)==0){
    next
  }
  tmpdfp8[ix,Novel_weight:= sum(p8$data$DistFromRoot.x[iy])]
  # print(ix)
}
p8$data$AnyKnown <- NULL
tmpdfp8$Weight_Novel = tmpdfp8$Novel_weight / tmpdfp8$Total_weight
tmpdfp8 <- distinct(tmpdfp8)
p8$data <- merge(p8$data,data.frame(tmpdfp8)[,c("node","Weight_Novel","AnyKnown")],by="node",all.x=T,all.y=F)
p8$data$AnyKnown[novel] <- F
p8$data$AnyKnown[old] <- T
p8$data$AnyKnown <- as.num(p8$data$AnyKnown)
w1 <- wh(p8$data$isTip)
p8$data$Weight_Novel[intersect(w1,wh(p8$data$AnyKnown==1))] = 0 
p8$data$Weight_Novel[setdiff(w1,wh(p8$data$AnyKnown==1))] = 1 #p8$data$DistFromRoot.x[setdiff(w1,wh(p8$data$AnyKnown==1))]


p9 <- p8 + geom_tree(mapping = aes(color = AnyKnown)) + scale_color_gradient(low="#565656", high="#22ace7")# +
# p9 <- p8 + scale_color_gradient(aes(color = AnyKnown),low="#565656", high="#22ace7")# +
QuickRDSSave(object = p9)
p9$data$Note[wh(p9$data$Note == "X")] = LysisSign
p9$data$Note[wh(p9$data$node %in% Putative_CRISPR_Matching_Nodes)] = CRISPRSign

p6 = p9 +   geom_tiplab(data = td_filter(!is.na(Note)), aes(image = Note),size = 0.005, geom = "image",align = T,hjust=.5,linesize=0.05) #+

Groups2Tag <- p6$data$Taxa[wh(p6$data$rank == "Order")]
Groups2Tag <- Groups2Tag[-grep(pattern = fixed("."),x = Groups2Tag,fixed = T)]
Groups2Tag <- unique(c(Groups2Tag,"Reverse_Transcriptases","Magsaviricetes","Matonaviridae","Hepelivirales","Pisuviricota","Duplornaviricota","p.0002.base-Kitrino","Negarnaviricota","f.0278","Lenarviricota","Kitrinoviricota","genPartiti.0019.base-Deltapartitivirus", "Cystoviridae","Picobirnaviridae","Partitiviridae","f.0271.base-Toga","f.0273.base-Toga","f.0268.base-Toga","f.0066.base-Hypo","f.0069.base-Hypo","f.0068.base-Hypo","f.0215.base-Deltaflexi","Tymoviridae","Nidovirales","Tobaniviridae","Hypoviridae","Deltaflexiviridae","Benyviridae","Closteroviridae","Botourmiaviridae","f.0024.base-Solinvi","Polycipiviridae","Rhabdoviridae","Coronaviridae"))
Groups2Tag <- setdiff(Groups2Tag,c("Muvirales"))
p6$data$NymN[wh(!(p6$data$Taxa %in% Groups2Tag))] = NA

p6$data$NymN <- gsub(x = p6$data$NymN,replacement = "",pattern = gsub(x = toString(RankDF$RegVec),pattern = ", ", replacement = "|"))
p6$data <- distinct(p6$data)
p11 <- p6 +  geom_label_repel(aes(label=NymN,fill=rank), segment.linetype = 6, segment.colour="Red", label.size = 0.01, max.time = 1, max.iter = 100, max.overlaps = 30, label.padding = 0.0, na.rm = T,
                              fontface="italic",
                              arrow = arrow(length = unit(0.01, "cm"),type="closed",ends = "last"))
# segment.curvature = -1e-5,
# min.segment.length = 1,
# ) #,force = 1,force_pull = 0,size=2.0,label.size = 0.01,nudge_x = 0.05,na.rm = T,segment.colour="RED",max.overlaps =1000)

### Putting labels on the Phyla
w1 <- filter(p11$data,rank == "Phylum")
w1 <- Rename1Col(w1,"Taxa","cladelab")
w1 <- Rename1Col(w1,"node","nodex")
w1$cladelab[wh(w1$cladelab=="p.0002.base-Kitrino")] <- "<----p.0002----|"
# colorcols2 <- colorcols
p12 <- p11 +#  scale_fill_manual(values = colorcols) + 
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Lenarviricota")],label =w1$cladelab[wh(w1$cladelab=="Lenarviricota")], barcolour=unname(colorcols["Lenarviricota"]), offset = -0.8, align = T, offset.text = .02,hjust = "center", barsize = 0.9, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Pisuviricota")], label =w1$cladelab[wh(w1$cladelab=="Pisuviricota")],barcolour=unname(colorcols["Pisuviricota"]), offset = -1.1, align = T, offset.text = .02,hjust = "center", barsize = 1.0, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Kitrinoviricota")], label =w1$cladelab[wh(w1$cladelab=="Kitrinoviricota")],barcolour=unname(colorcols["Kitrinoviricota"]), align = T,offset = -0.5, offset.text = .02,hjust = "center", barsize = 1.0, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Duplornaviricota")], label =w1$cladelab[wh(w1$cladelab=="Duplornaviricota")],barcolour=unname(colorcols["Duplornaviricota"]), align = T, offset = -1, offset.text = .02,hjust = "center", barsize = .5, fontsize = 4, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Negarnaviricota")], label =w1$cladelab[wh(w1$cladelab=="Negarnaviricota")],barcolour=unname(colorcols["Negarnaviricota"]), align = T, offset = -1.2, offset.text = .02,hjust = "center", barsize = .8, fontsize = 6, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="<----p.0002----|")], label =w1$cladelab[wh(w1$cladelab=="<----p.0002----|")],barcolour=unname(colorcols["p.0002"]), align = T, offset = -1, offset.text = .1,hjust = "center", barsize = .5, fontsize = 6, angle = 330, horizontal = F, fontface="italic") 


p12$data$Alt=NA
p12$data$Alt[wh(p12$data$node %in% AltSets$node)] <- T

p12$data$Perm=NA
p12$data$Perm[wh(p12$data$node %in% PermutSets$node)] <- T

p12$data$shapo = NA
p12$data$shapo[wh(p12$data$Alt)] <- "Alternative Genetic Code"
p12$data$shapo[wh(p12$data$Perm)] <- "Permuted RdRP"

p12$data$NkidsProp[] <-(16+3*(scales::rescale(x = as.integer(p12$data$Nkids),to=c(1,2))))

p13 <- p12 + geom_point2(aes(size=(6+(NkidsProp*24)),shape=shapo), na.rm = T,color = "purple", fill='purple')
p13 <- p12 +   geom_nodepoint(mapping = aes(subset = (!is.na(shapo)),shape=shapo,size = (8+(NkidsProp*24)),fill=shapo),)#,color = "purple", fill='purple') #+

# QuickSave
tmppppp <- p13
tmppppp$layers[[1]] = NULL
tmppppp$data <- `dropcols<-`(tmppppp$data,c("DistFromRoot.x","Lysis.y","Lysis.x","CRISPR.x","CRISPR.y","DistFromRoot.y"))
fastSave::saveRDS.lbzip2(n.cores = 12,object = tmppppp,file = p0("AlignedTip.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".RDS"))
{
  source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
  setwd("/media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/workdir/")
  library(ggtreeExtra);library(ggstar);library(ggnewscale);library(ggtree);library(stringr)
  rm(testasd)
  testasd <- QuickRDSRead(n.cores=11,fname="tmp.*.RDS")
  testasd <- fastSave::readRDS.lbzip2(n.cores = 12,file = list.files(path = ".",pattern = "^AlignedTip.*RDS"))
  pdf(encoding = "ISOLatin1.enc",file = paste0("AlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.pdf"),width = 41, height = 41,title = "Riboviria (CRISPR+Lysis stared, Permutation + Alt code as internal shapes)",useDingbats = T)
  testasd# +theme(legend.position="bottom")
  dev.off()
  
}
agg_png("tmptest.png", width = 1000, height = 500, res = 144,)
Cairo::CairoSVG(antialias = "none",file = paste0("SVGAlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.svg"),width = 45, height = 45)
Cairo::CairoSVG(antialias = "none",file = paste0("SVGAlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.svg"),width = 45, height = 45)
testasd# +theme(legend.position="bottom")
dev.off()
Cairo::CairoSVG(antialias = "none",file = paste0("cSVG.AlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.svg"),width = 45, height = 45)
Cairo::CairoSVG(file = paste0("cSVG.AlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.svg"),width = 45, height = 45,)
Cairo::CairoPNG(filename = paste0("cPNG.AlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.svg"),width = 45, height = 45)
CairoPNG(filename = "Rplot%03d.png", width = 480, height = 480,
         pointsize = 12, bg = "white",  res = NA, ...)


QuickTreePDF()
tmppppp#+theme(legend.position="bottom")
dev.off()

# agg_png("tmptest.png", width = 1000, height = 500, res = 144)
# cairo_pdf(antialias = "none")
pdf(file = paste0("p71.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.pdf"),width = 41, height = 41,title = "Riboviria (CRISPR+Lysis stared)",useDingbats = T)

# 
# tt = '((snail,mushroom),(((sunflower,evergreen_tree),leaves),green_salad));'
# tree = read.tree(text = tt)
# d <- data.frame(label = c('snail','mushroom', 'sunflower',
#                           'evergreen_tree','leaves', 'green_salad'),
#                 group = c('animal', 'fungi', 'flowering plant',
#                           'conifers', 'ferns', 'mosses'))
# # 
# ggtree(tree, linetype = "dashed", size=1, color='firebrick') %<+% d + 
#   xlim(0, 4.5) + ylim(0.5, 6.5) +
#   geom_tiplab(parse="emoji", size=15, vjust=.25) +
#   geom_tiplab(aes(label = group), geom="label", x=3.5, hjust=.5)
# 
# geom_cladelabel(node=53, label="duck", parse="emoji",
#                 fontsize=12, align=TRUE, colour="steelblue") +
# p6 <-   p5 + geom_cladelabel(node=97141, label="duck", parse="emoji",
#                   fontsize=12, align=TRUE, colour="darkkhaki")
# p5 + geom_cladelabel(node=41, label="bacteria", parse="emoji",
#                 fontsize=12, align=TRUE, colour="firebrick")
# 
# 
# tr <- rtree(30)
# 
# ggtree(tr) + xlim(NA, 5.2)  +
#   geom_cladelabel(node=53, label="phage", parse="emoji",
#                   fontsize=12, align=TRUE, colour="steelblue") +
#   geom_cladelabel(node=48, label="family", parse="emoji",
#                   fontsize=12, align=TRUE, colour="darkkhaki")
# emoji('dog')
# 
# 

# info2 = filter(TaxDT4W, Rank %in% c("Class","Order","Family"))
# info2 = filter(info2, MaxQtc > 0.5)
# info2  = Rename1Col(info2,"MaxQtc_Node","node")
# info2 = info2[, c("Taxa"  , "node" , "Rank")]
# p5 = p3
# p5$data <- merge((p5$data),
#                  info2,
#                  by = "node",
#                  all.x = T,
#                  all.y = F)

# p6 = p5 +  geom_label_repel(aes(label=Taxa),force = 1,force_pull = 0,size=2.0,label.size = 0.01,nudge_x = 0.05,na.rm = T,segment.colour="RED",max.overlaps =1000)

# pdf(
#   "TmpTree.RAW.pdf",
#   width = 20,
#   height = 20,
#   title = "RAW"
# )
QuickTreePDF()
# p4
# p3
p6
# p7
dev.off()

###### Reduced tree ###### 
tmptree <- WorkTree
infoInternals$Nkids <- lengths((GetKids(WorkTree,as.num(infoInternals$node),type="tips")))
infoInternals$NkidsProp <- infoInternals$Nkids/Ntip(WorkTree)
Fams2Plot <- data.frame("Name"=infoInternals$Taxa, BestKid="") 
Fams2Plot$BestKids = vector('list',c(nrow(Fams2Plot)))   # List place holder 1
w1 <- wh(Fams2Plot$Name %in% c(infoInternals$Taxa[wh(infoInternals$rank=="Genus")]))
w1 <- setdiff(w1,wh(Fams2Plot$Name == "genPartiti.0019.base-Deltapartitivirus"))
Fams2Plot <- Fams2Plot[-w1,]
# Fams2Plot <- setdiff(c(infoInternals$Taxa[infoInternals$rank == "Family"],"Reverse_Transcriptases","genPartiti.0019.base-Deltapartitivirus"),c("Blumeviridae","Alvernaviridae", "Barnaviridae", "Carmotetraviridae", "Quadriviridae",  "Sunviridae"))
TaxaScaler <- 3
for(iy in 1:nrow(Fams2Plot)){
  ix <- Fams2Plot$Name[iy]
  FamNodeLabel <- (infoInternals$node.label[wh(infoInternals$Taxa == ix)])
  FamNodeid <- NLabel2NID(tmptree,FamNodeLabel)
  if(!(FamNodeLabel %in% c(tmptree$node.label,tmptree$tip.label))){
    next
  }
  # Pick the tip that's furthest from the root.
  tmpdf <- filter(info2,node.label %in% unlist(GetKids(tmptree,FamNodeid,Return_Labels = T)))
  setorderv(setDT(tmpdf),c("DistFromRoot"), order = c(-1))
  BestKid <- tmpdf$node.label[1]
  Fams2Plot$BestKid[iy] <- BestKid
  Ngets <- (TaxaScaler*(1+infoInternals$NkidsProp[wh(infoInternals$Taxa == ix)]))
  Fams2Plot$BestKids[iy] <- list(sample(x = tmpdf$node.label,size = round(Ngets),replace = T))
}
# wt1 <- p11$data$label[intersect(wh(p11$data$Note==LysisSign),wh(p11$data$Phylum!="Lenarviricota")) ]
# wt2 <- p11$data$label[wh(p11$data$Note==CRISPRSign) ]

tips2keep <- c(wt2,wt1,"rv20_v2_2132","Rv4_075436","rv20_v2_5848","Rv4_319121","Rv4_314863",unique(unlist(unlist(Fams2Plot$BestKids))))
tmptree <- drop.tip(phy = WorkTree,tip = setdiff(WorkTree$tip.label,unique(c(tips2keep,Fams2Plot$BestKid))),trim.internal = T,collapse.singles = F)

tmpestdf <- filter(p12$data,label %in% reduced_p4$data$label)[,c("label","NkidsProp","isTip")]
# tmpestdf <- tmpestdf[!tmpestdf$isTip,]
tmpest4map <- (tmpestdf$NkidsProp+1) * 360
names(tmpest4map) <- tmpestdf$label

reduced_p4 <- ggtree(tmptree,layout = "circular", mapping = aes(x ~ x*tmpestdf$NkidsProp))
reduced_p4 <- scaleClade(ggtree(tmptree,layout = "circular"),node = 6625,scale = 2.81,vertical_only = T) # Lena
reduced_p4 <- scaleClade(reduced_p4,node = 1389,scale = 0.945,vertical_only = T) # Pisu
reduced_p4 <- scaleClade(reduced_p4,node = 3471,scale = 0.80,vertical_only = T) # Kitrino
reduced_p4 <- scaleClade(reduced_p4,node = 5642,scale = 1.97,vertical_only =T) # Duplo
reduced_p4 <- scaleClade(reduced_p4,node = 5939,scale = 0.20,vertical_only = T) # Nega
reduced_p4 <- scaleClade(reduced_p4,node = 7726,scale = 0.059,vertical_only =T) # RTs
reduced_p4 <- scaleClade(reduced_p4,node = 5618,scale = 0.02,vertical_only = T) # p.0001
reduced_p4 <- scaleClade(reduced_p4,node = 5620,scale = 0.06,vertical_only =T) # p.0002


reduced_p12$data[wh(reduced_p12$data$rank == "Phylum"),c("node","NkidsProp","Taxa")]
p12$data[wh(p12$data$rank == "Phylum"),c("node","NkidsProp","Taxa")]

reduced_p12$data[wh(reduced_p12$data$Taxa == "Kitrinoviricota"),c("node","NkidsProp")]


reduced_p4$data <- merge((reduced_p4$data),
                         Rename1Col(`dropcols<-`(infoInternals,"node"),"node.label","label"),
                         by = "label",
                         all.x = T,
                         all.y = F)


reduced_p4$data <- merge((reduced_p4$data),
                         Rename1Col(`dropcols<-`(info2,c("node","Lysis","CRISPR","DistFromRoot")),"node.label","label"),
                         by = "label",
                         all.x = T,
                         all.y = F)


tmpinfoInternals <- infoInternals
tmpinfoInternals$AllKids <-  vector('list',c(nrow(tmpinfoInternals)))
tmpinfoInternals$newnode <- 0
for(i in 1:nrow(tmpinfoInternals)){
  tmpid <- NLabel2NID(tmptree,tmpinfoInternals$node.label[i])
  if(is.na(tmpid)){
    print(p0(tmpinfoInternals$Taxa[i]," Isn't on the tree"))
    next
  }
  tmpinfoInternals$newnode[i] <- tmpid
  tmpkids <- unlist(GetKids(tmptree,tmpid,Return_Labels = T,type="all"))
  tmpinfoInternals$AllKids[i] <- list(tmpkids)
  w1 <- wh(reduced_p4$data$label %in% tmpkids)
  reduced_p4$data[w1,tmpinfoInternals$rank[i]] <- tmpinfoInternals$Taxa[i]
}



reduced_p7 <- reduced_p4 +scale_fill_manual(values = colorcols)+
  geom_fruit(geom=geom_bar,mapping=aes(x=Mlen,fill=Phylum),
             pwidth=0.5,orientation="y",stat="identity",
             grid.params=list(vline =F),
             axis.params=list(title = "Maximal observed Genome Length [nt]",
                              axis       = "x",
                              title.height = 0.001,
                              text.size = 3,
                              hjust      = -0.2,
                              vjust      = +0.1,
                              nbreak     = 10,
                              line.size = 0.4)) +
  theme(axis.title = element_text(hjust = -5,vjust = -4))
reduced_p7[["layers"]][[7]][["data"]][["y"]] <- 4
# 
# AltSets <- fread("Alt_Clades.tab",sep = "\t")
# AltSets$leaflist = vector('list',c(nrow(AltSets)))   # List place holder 1
# AltSets$leaflist <- apply(AltSets,MARGIN = 1,FUN = function(x) as.char(unlist(str_split_fixed(string = x["V2"],pattern = ",",n = Inf))))
# AltSets$node <- unlist(apply(AltSets,MARGIN = 1,FUN = function(x) get_mrca_of_set(WorkTree,unlist(x["leaflist"]))))
# AltSets$node.label <- NID2NLabel(WorkTree,AltSets$node)
# 
# PermutSets <- fread("Perm_Clades.txt",sep = "\t")
# PermutSets$leaflist = vector('list',c(nrow(PermutSets)))   # List place holder 1
# PermutSets$leaflist <- apply(PermutSets,MARGIN = 1,FUN = function(x) as.char(unlist(str_split_fixed(string = x["V2"],pattern = ",",n = Inf))))
# PermutSets$node <- unlist(apply(PermutSets,MARGIN = 1,FUN = function(x) get_mrca_of_set(WorkTree,unlist(x["leaflist"]))))
# PermutSets$node.label <- NID2NLabel(WorkTree,PermutSets$node)
# intersect(PermutSets$node,AltSets$node)
# AltSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/Alt_Sign.png")
# PermuteSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/Permute_Sign.png")
# CRISPRSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/BlueStar.png")
# LysisSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/RedStar.png")

reduced_p71  <- reduced_p7
reduced_p71$data$Nkids <- lengths((GetKids(tmptree,reduced_p71$data$node,type="tips")))
reduced_p71$data$NkidsProp <- reduced_p71$data$Nkids/Ntip(tmptree)
reduced_p71$data$shapo = NA
reduced_p71$data$shapo[wh(reduced_p71$data$label %in% AltSets$node.label)] <- "17"
reduced_p71$data$shapo[wh(reduced_p71$data$label %in% PermutSets$node.label)] <- "15"
reduced_p71 <- reduced_p71 + geom_point2(aes(subset = (label %in% c(AltSets$node.label,PermutSets$node.label)),size=(6+(NkidsProp*28)),shape=shapo), na.rm = T,color = "purple", fill='purple')
# reduced_p71 <- reduced_p71 + geom_point2(aes(subset = (label %in% PermutSets$node.label),size=6+(NkidsProp*28)),na.rm = T, shape=15,  color = "purple", fill='purple')# ,size=0.5)
# fastSave::saveRDS.lbzip2(n.cores = 12,object = reduced_p71,file = p0("reduced_p71",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".RDS"))
# pdf(width = 41,height = 41)
# reduced_p71
# dev.off()


reduced_p8 <- reduced_p71

# 
# Nodes2Iter <- wna(p8$data$Novel)
# novel <- wh(p8$data$Novel)
# old <- wh(!(p8$data$Novel))
# p8$data$Novel_weight <- 0
# p8$data$Novel_weight[novel] <- 1
# p8$data$AnyKnown <- T
# p8$data$AnyKnown <- !p8$data$Novel
# tmpdfp8 <- data.frame(p8$data[Nodes2Iter,c("node","label","Novel_weight","DistFromRoot.x","AnyKnown")])
# tmpdfp8$Total_weight <- 0
# tmpdfp8$AllKids = GetKids(WorkTree,tmpdfp8$node,type="tips")
# setDT(tmpdfp8)
# for(ix in 1:nrow(tmpdfp8)){
#   iy <- wh(p8$data$node %in% unlist(tmpdfp8$AllKids[ix]))
#   tmpdfp8[ix,Total_weight:= sum(p8$data$DistFromRoot.x[iy])]
#   tmpdfp8[ix,AnyKnown:= AnyFalse(p8$data$Novel[iy])]
#   iy <- wh(p8$data$Novel[iy] ==1)
#   if(len(iy)==0){
#     next
#   }
#   tmpdfp8[ix,Novel_weight:= sum(p8$data$DistFromRoot.x[iy])]
#     # print(ix)
# }
# p8$data$AnyKnown <- NULL
# tmpdfp8$Weight_Novel = tmpdfp8$Novel_weight / tmpdfp8$Total_weight

reduced_p8$data <- merge(reduced_p8$data,data.frame(tmpdfp8)[,c("label","Weight_Novel","AnyKnown")],by="label",all.x=T,all.y=F)
reduced_p8$data$AnyKnown[wh(reduced_p8$data$Novel)] <- F
reduced_p8$data$AnyKnown[ wh(!(reduced_p8$data$Novel))] <- T
reduced_p8$data$AnyKnown <- as.num(reduced_p8$data$AnyKnown)
w1 <- wh(reduced_p8$data$isTip)
reduced_p8$data$Weight_Novel[intersect(w1,wh(reduced_p8$data$AnyKnown==1))] = 0 
reduced_p8$data$Weight_Novel[setdiff(w1,wh(reduced_p8$data$AnyKnown==1))] = 1 #p8$data$DistFromRoot.x[setdiff(w1,wh(p8$data$AnyKnown==1))]


reduced_p9 <- reduced_p8 + geom_tree(mapping = aes(color = AnyKnown)) + scale_color_gradient(low="#565656", high="#19A2DD")# +
reduced_p9$data$Note[wh(reduced_p9$data$Note == "X")] = LysisSign
reduced_p9$data$Note[wh(reduced_p9$data$node %in% Putative_CRISPR_Matching_Nodes)] = CRISPRSign

reduced_p6 = reduced_p9 +   geom_tiplab(data = td_filter(!is.na(Note)), aes(image = Note),size = 0.005, geom = "image",align = T,hjust=.5,linesize=0.05) #+

Groups2Tag <- reduced_p6$data$Taxa[wh(reduced_p6$data$rank == "Order")]
Groups2Tag <- Groups2Tag[-grep(pattern = fixed("."),x = Groups2Tag,fixed = T)]
Groups2Tag <- unique(c(Groups2Tag,"Reverse_Transcriptases","Magsaviricetes","Matonaviridae","Hepelivirales","Pisuviricota","Duplornaviricota","p.0002.base-Kitrino","Negarnaviricota","f.0278","Lenarviricota","Kitrinoviricota","genPartiti.0019.base-Deltapartitivirus", "Cystoviridae","Picobirnaviridae","Partitiviridae","f.0271.base-Toga","f.0273.base-Toga","f.0268.base-Toga","f.0066.base-Hypo","f.0069.base-Hypo","f.0068.base-Hypo","f.0215.base-Deltaflexi","Tymoviridae","Nidovirales","Tobaniviridae","Hypoviridae","Deltaflexiviridae","Benyviridae","Closteroviridae","Botourmiaviridae","f.0024.base-Solinvi","Polycipiviridae","Rhabdoviridae","Coronaviridae"))
Groups2Tag <- setdiff(Groups2Tag,c("Muvirales"))
reduced_p6$data$NymN[wh(!(reduced_p6$data$Taxa %in% Groups2Tag))] = NA

reduced_p6$data$NymN <- gsub(x = reduced_p6$data$NymN,replacement = "",pattern = gsub(x = toString(RankDF$RegVec),pattern = ", ", replacement = "|"))
reduced_p6$data <- distinct(reduced_p6$data)
reduced_p11 <- reduced_p6 +  geom_label_repel(aes(label=NymN,fill=rank), segment.linetype = 6, segment.colour="Red", label.size = 0.01, max.time = 1, max.iter = 100, max.overlaps = 30, label.padding = 0.0, na.rm = T,
                                              fontface="italic",
                                              arrow = arrow(length = unit(0.01, "cm"),type="closed",ends = "last"))

### Putting labels on the Phyla
w1 <- filter(reduced_p11$data,rank == "Phylum")
w1 <- Rename1Col(w1,"Taxa","cladelab")
w1 <- Rename1Col(w1,"node","nodex")
w1$cladelab[wh(w1$cladelab=="p.0002.base-Kitrino")] <- "p.0002"
reduced_p12 <- reduced_p11 +#  scale_fill_manual(values = colorcols) + 
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Lenarviricota")],label =w1$cladelab[wh(w1$cladelab=="Lenarviricota")], barcolour=unname(colorcols["Lenarviricota"]), offset = -0.8, align = T, offset.text = .001,hjust = "center", barsize = 2.4, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Pisuviricota")], label =w1$cladelab[wh(w1$cladelab=="Pisuviricota")],barcolour=unname(colorcols["Pisuviricota"]), offset = -0.9, align = T, offset.text = .001,hjust = "center", barsize = 2.4, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Kitrinoviricota")], label =w1$cladelab[wh(w1$cladelab=="Kitrinoviricota")],barcolour=unname(colorcols["Kitrinoviricota"]), align = T,offset = -0.5, offset.text = .001,hjust = "center", barsize = 2.4, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Duplornaviricota")], label =w1$cladelab[wh(w1$cladelab=="Duplornaviricota")],barcolour=unname(colorcols["Duplornaviricota"]), align = T, offset = -1, offset.text = .001,hjust = "center", barsize = 2.4, fontsize = 4, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Negarnaviricota")], label =w1$cladelab[wh(w1$cladelab=="Negarnaviricota")],barcolour=unname(colorcols["Negarnaviricota"]), align = T, offset = -1.2, offset.text = .001,hjust = "center", barsize = 2.4, fontsize = 6, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="p.0002")], label =w1$cladelab[wh(w1$cladelab=="p.0002")],barcolour=unname(colorcols["p.0002"]), align = T, offset = -1, offset.text = .1,hjust = "center", barsize = 3.2, fontsize = 6, angle = 30, horizontal = F, fontface="italic") 


# reduced_p12$data <- distinct(`dropcols<-`(reduced_p12$data)
QuickTreePDF()
scaleClade(reduced_p12,node = reduced_p12$data$node[wh(reduced_p12$data$Taxa == "Lenarviricota")],scale = 0.3)
dev.off()
# QuickSave
tmppppp <- reduced_p12
tmppppp$layers[[1]] = NULL
tmppppp$data <- `dropcols<-`(tmppppp$data,c("DistFromRoot.x","Lysis.y","Lysis.x","CRISPR.x","CRISPR.y","DistFromRoot.y"))
fastSave::saveRDS.lbzip2(n.cores = 12,object = tmppppp,file = p0("AlignedTip.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".RDS"))
{
  setwd("/media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/workdir/")
  library(ggtreeExtra);library(ggstar);library(ggnewscale);library(ggtree);library(stringr)
  rm(testasd)
  testasd <- fastSave::readRDS.lbzip2(n.cores = 12,file = list.files(path = ".",pattern = "^AlignedTip.*RDS"))
  pdf(file = paste0("AlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.pdf"),width = 41, height = 41,title = "Riboviria (CRISPR+Lysis stared)",useDingbats = T)
  testasd# +theme(legend.position="bottom")
  dev.off()
  
}




p3 = p4 
p3 = p3 +   geom_tiplab(aes(label = Note), geom="label", x=3.5, hjust=.5,align = F) +
  ggtitle(label = "'X's are putative lysis encoding, 'C's are putative CRISPR matching.") +
  theme(plot.title = element_text(size = 50, hjust = 0.5,face = "bold"))

reduced_p6 = p3
View(reduced_p8$data)
Groups2Tag <- reduced_p6$data$Taxa[wh(reduced_p6$data$rank == "Order")]
Groups2Tag <- c(Groups2Tag,"Reverse_Transcriptases","genPartiti.0019.base-Deltapartitivirus", "Cystoviridae","Picobirnaviridae","Partitiviridae")

reduced_p6$data$NymN[wh(!(reduced_p6$data$Taxa %in% Groups2Tag))] = NA
reduced_p6 = reduced_p6 +  geom_label_repel(aes(label=NymN)) #,force = 1,force_pull = 0,size=2.0,label.size = 0.01,nudge_x = 0.05,na.rm = T,segment.colour="RED",max.overlaps =1000)

reduced_p8 <- reduced_p6

# Nodes2Iter <- len(wh(!(p8$data$isTip)))
Nodes2Iter <- reduced_p8$data$label[wna(reduced_p8$data$Novel)]
reduced_p8$data$Novel[wh(reduced_p8$data$Novel)] <- 1
reduced_p8$data$Novel[wh(!(reduced_p8$data$Novel))] <- 0
tmpdf <- data.frame(label=Nodes2Iter,Novel=0)

tmpdf$AllKids = GetKids(tmptree,NLabel2NID(tmptree,Nodes2Iter))
for(ix in 1:nrow(tmpdf)){
  iy <- wh(reduced_p8$data$node %in% unlist(tmpdf$AllKids[ix]))
  tmpdf$Novel[ix] <- mean(reduced_p8$data$Novel[iy])
  reduced_p8$data$Novel[wh(reduced_p8$data$label == tmpdf$label[ix])] <-  mean(reduced_p8$data$Novel[iy])
  
}

reduced_p9 <- reduced_p8 + geom_tree(mapping = aes(color = Novel)) +scale_color_continuous(low="black", high="cyan")# +

# View(p9$data)

p7 <- reduced_p9  + 
  # new_scale_fill() + 
  geom_fruit(geom=geom_bar,mapping=aes(x=Mlen,fill=Class),
             pwidth=0.5,orientation="y",stat="identity",
             grid.params=list(vline =F,),
             axis.params=list(title = "Maximal observed Genome Length ",
                              axis       = "x",
                              title.height = 0.001,
                              text.size = 2,
                              hjust      = 1,
                              vjust      = 0.5,
                              nbreak     = 6,
             ))

QuickTreePDF()
p7+theme(legend.position="bottom") 
dev.off()






