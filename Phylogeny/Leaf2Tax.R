# Information  --------------------------------------------------------------------
## Script name: Leaf2Tax.R
## Description: Affiliate tree leaves (RCR90) taxonomic labels (ICTV MSL / NCBI's Taxonomy).
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

#### 17.10.2021 ####
AllSneakies = unique(scan("Misc/All_Sneakies.lst", what = "\n"))
AllSneakies2 = stringi::stri_split_fixed(pattern = "|",  n = 5,  str = AllSneakies,  simplify = T)[, 4]

VMR_MSL36 = ReadXlsx("/media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/Misc/VMR 200721 MSL36.xlsx")
VMR_MSL36 = filter(VMR_MSL36, Realm == "Riboviria")# !(Realm %in% c("Varidnaviria","Adnaviria","Monodnaviria", "Duplodnaviria")))
VMR_MSL36 = RemoveNACols(VMR_MSL36)
VMR_MSL36 <- `dropcols<-`(VMR_MSL36,c(1,2))
colnames(VMR_MSL36) = gsub(pattern = " ",replacement = "_",x = colnames(VMR_MSL36))
VMR_MSL36$RefSeqs = vector('list',c(nrow(VMR_MSL36)))   # List place holder 1
VMR_MSL36$RefSeqs = stringr::str_extract_all(string = (VMR_MSL36$Virus_REFSEQ_accession),pattern =  "NC_([0-9]+)")
VMR_MSL36 = VMR_MSL36[-wna(VMR_MSL36$RefSeqs),]

writeLines(unlist(VMR_MSL36$RefSeqs),"RefSeqs2Get.txt")

RefSeqs2Get = readDNAStringSet("RefSeqs2Get.fasta")
RefSeqs2GetDF = XString2DF(faa = RefSeqs2Get,input_was = "RefSeq",trimwhite = T,seqcolname = "seq",addlength = T)
names(RefSeqs2Get) = str_split_fixed(string = names(RefSeqs2Get), pattern = " ", n = Inf)[,1]
RefSeqs2GetDF$RefSeqID.version = names(RefSeqs2Get) 
RefSeqs2GetDF$RefSeqID = tstrsplit(x =RefSeqs2GetDF$RefSeqID,split = ".",fixed =T )[[1]]
names(RefSeqs2Get) =  RefSeqs2GetDF$RefSeqID 
writeXStringSet(RefSeqs2Get,"RefSeqs2Get.fas")

VMR_MSL36$RefSeqsChar = apply(X = VMR_MSL36,MARGIN = 1, FUN = function(x) toString(unlist(unname(x["RefSeqs"]))))
VMR_MSL36$RefSeqs = NULL
VMR_MSL36 = TabUnroll(dfdt = VMR_MSL36,sep = ", ",colnym = "RefSeqsChar","RefSeqID")

Blastn_IDFT_vs_ICTV = fread("../Valerian/IDFT_blastn_ICTV/blastn_IDFT_vs_ICTV_out.tsv", col.names =c("ND", "RefSeqID", "qcov", "nid", "pos", "pid", "qL", "pL",  "q1", "q2", "p1", "p2", "length", "evalue", "score", "gaps"))
DCB <- Blastn_IDFT_vs_ICTV 
DCB = CalcPcoverage(DCB)
DCB = CalcQcoverage(DCB)

DCB <- filter(DCB, ND %in% IDDT5$ND)

w1 = unique(intersect(which(DCB$qcov >= 95),which(DCB$nid >= 100)))
DCB = DCB[w1,]
DCB = filter(DCB, pid >= 98 )
DCB = filter(DCB, evalue <= 1e-30)
DCB = filter(DCB, length > 200)
DCB = filter(DCB, qCoverage > 0.75)
DCB = filter(DCB, length/qL > 0.95)

setorderv(setDT(DCB),c("ND", "score", "nid", "qcov", "evalue"), order = c(-1,-1,-1,-1,1))
DCB = unique(DCB, incomparables = FALSE, fromLast = FALSE, by="ND")
DCB = merge(DCB,`dropcols<-`(VMR_MSL36,c("RefSeqsChar","Genome_coverage")),by="RefSeqID", all.x=T, all.y=F)
DCB = merge(DCB,IDDT5[,c("ND","RID")],by="ND",all.x=T,all.y=F)


TB.lst = scan("Misc/TB.lst",what = "\n")
NB.lst = setdiff(scan("Misc/NB.lst",what = "\n"),TB.lst)
Rnkcls = c("Class","Kingdom","Order","Family","Genus","Species","Phylum")

nr210930.id2tt.tab = fread("/media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/nr210930.id2tt.tab", col.names = c("RID", "GbID", "taxid", "lin"))
taxiddf_lin = data.frame(taxonomizr::getTaxonomy(ids = nr210930.id2tt.tab$taxid, sqlFile = "/media/HDD1/uri/DBs/NCBI_Tax/accessionTaxa.sql",
                                                 desiredTaxa = c("subkingdom",
                                                                 "subclass",
                                                                 "suborder",
                                                                 "subfamily",
                                                                 "subphylum",
                                                                 "species",
                                                                 "genus",
                                                                 "no rank",
                                                                 "family",
                                                                 "order",
                                                                 "class",
                                                                 "phylum",
                                                                 "kingdom",
                                                                 "clade",
                                                                 "superkingdom",
                                                                 "subgenus")))
taxiddf_lin = RemovesSparseols(taxiddf_lin)
taxiddf_lin = bind_cols(nr210930.id2tt.tab, taxiddf_lin)
WriteWolfTbl(taxiddf_lin, "taxiddf_lin.tab")
taxiddf_lin$class[which(taxiddf_lin$class == "Allassoviricetes")] = "Leviviricetes"

Levi2021prop = fread0("Levi2021prop.tsv")
colnames(Levi2021prop) = tolower(colnames(Levi2021prop))
Levi2021prop$taxid = taxonomizr::accessionToTaxa(p0(Levi2021prop$gbid, ".1"),
                                                 '/media/HDD1/uri/DBs/NCBI_Tax/accessionTaxa.sql')
Levi2021prop = filter(Levi2021prop, (taxid %in% taxiddf_lin$taxid))
Levi2021prop = merge(distinct(data.frame(Levi2021prop)), taxiddf_lin[, c("taxid", "RID", "GbID")], by = "taxid", all.x = T, all.y = F)
taxiddf_lin1 = filter(taxiddf_lin, !(taxid %in% Levi2021prop$taxid))
taxiddf_lin1 = distinct(plyr::rbind.fill(taxiddf_lin1, Levi2021prop))

#### add_lev -->Lenarviricota, add_pla --> Pisuviricota;Stelpaviricetes;Stellavirales ####
taxiddf_lin1$superkingdom[grep(x = taxiddf_lin1$RID, pattern = "al2_", fixed = T)] = "Viruses"
taxiddf_lin1$superkingdom[grep(x = taxiddf_lin1$RID, pattern = "ap2_", fixed = T)] = "Viruses"
taxiddf_lin1$phylum[grep(x = taxiddf_lin1$RID, pattern = "al2_", fixed = T)] = "Lenarviricota"
taxiddf_lin1$phylum[grep(x = taxiddf_lin1$RID, pattern = "ap2_", fixed = T)] = "Pisuviricota"
taxiddf_lin1$class[grep(x = taxiddf_lin1$RID, pattern = "ap2_", fixed = T)] = "Stelpaviricetes"
taxiddf_lin1$order[grep(x = taxiddf_lin1$RID, pattern = "ap2_", fixed = T)] = "Stellavirales"

for (i in 1:ncol(taxiddf_lin1)){
  w1 = which(taxiddf_lin1[,i] == "")
  taxiddf_lin1[w1,i] = NA
}

for (i in c("phylum","order","family","class")){
  taxiddf_lin1[,i] = gsub(pattern = "unclassified ",ignore.case = T,replacement = "",x =taxiddf_lin1[,i] )
}

taxiddf_lin1$Nna = 0
taxiddf_lin1$Nna = apply(taxiddf_lin1[,c(5:19)], MARGIN = 1, FUN = function(x) (length(wna(x[]))))

tmpmat = str_split_fixed(string = taxiddf_lin1$lin,pattern = ";", n=Inf)

taxiddf_lin1[,RankDF$Rank] = NA
setDF(taxiddf_lin1)
for (i in 1:nrow(RankDF)){
  rnk = RankDF$Rank[i]
  rgx = RankDF$RegVec[i]
  w1  = (unlist(apply(tmpmat,MARGIN = 2, FUN = function(y) unlist(grep(pattern = rgx, x = y,value = F)))))
  if (class(w1) == "matrix"){
    print(w1)
    next
  }
  w2 = (unlist(apply(tmpmat,MARGIN = 2, FUN = function(y) unlist(grep(pattern = rgx, x = y,value = T)))))
  taxiddf_lin1[w1,rnk] = w2
}
for (i in c("Phylum","Order","Family","Class")){
  taxiddf_lin1[,i] = gsub(pattern = "unclassified ",ignore.case = T,replacement = "",x =taxiddf_lin1[,i] )
}

taxiddf_lin1[,toupper(RankDF$Rank)] = F
taxiddf_lin1[,setdiff(tolower(RankDF$Rank),colnames(taxiddf_lin1))] = NA
taxiddf_lin1[,p0("Hybrid_",toupper(RankDF$Rank))] = taxiddf_lin1[,tolower(RankDF$Rank)]

for (i in 1:nrow(RankDF)){
  rn2cmp = c(RankDF$Rank[i], tolower(RankDF$Rank[i]),toupper(RankDF$Rank[i])) 
  c1 = which(colnames(taxiddf_lin1) == rn2cmp[2] )
  if (length(c1) == 0 ){
    print(c1)
    next
  }
  w1  = which(taxiddf_lin1[,rn2cmp[1]] == taxiddf_lin1[,rn2cmp[2]])
  w1 = unique(w1,intersect(wna(taxiddf_lin1[,rn2cmp[1]]),wna(taxiddf_lin1[,rn2cmp[2]])))
  print(w1)
  
  if (length(w1) == 0 ){
    next
  }
  taxiddf_lin1[w1,rn2cmp[3]] = T
  
  taxiddf_lin1[w1,p0("Hybrid_",rn2cmp[3])] = taxiddf_lin1[w1,rn2cmp[1]] 
  w3 = setdiff(wna(taxiddf_lin1[,rn2cmp[1]]),wna(taxiddf_lin1[,rn2cmp[2]]))
  taxiddf_lin1[w3,p0("Hybrid_",rn2cmp[3])] = taxiddf_lin1[w3,rn2cmp[2]] 
  w4 = setdiff(wna(taxiddf_lin1[,rn2cmp[2]]),wna(taxiddf_lin1[,rn2cmp[1]]))
  taxiddf_lin1[w4,p0("Hybrid_",rn2cmp[3])] = taxiddf_lin1[w4,rn2cmp[1]] 
}

# Inspect to see if any sub-group label occurred in more than one parent group:
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Hybrid_CLASS","Hybrid_PHYLUM")]),vars ="Hybrid_CLASS" )
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Hybrid_CLASS","Hybrid_ORDER")]),vars ="Hybrid_ORDER" )
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Hybrid_FAMILY","Hybrid_ORDER")]),vars ="Hybrid_FAMILY" )
tmpdf = plyr::count(distinct(taxiddf_lin1[,c("Hybrid_FAMILY","Hybrid_GENUS")]),vars ="Hybrid_GENUS" )

TDF <- taxiddf_lin1[,unique(c("GbID","RID","species",p0("Hybrid_",toupper(RankDF$Rank))))]
TDF = RemoveNACols(TDF)
colnames(TDF) = gsub(pattern = "Hybrid_",ignore.case = T,replacement = "",x = colnames(TDF))
colnames(TDF) = tolower(colnames(TDF))
TDF$realm = "Riboviria"
TDF$superkingdom = "Viruses"
TDF$kingdom = "Orthornavirae"
TDF = distinct(TDF)
TDF$Nna = 0
TDF$Nna = apply(TDF, MARGIN = 1, FUN = function(x) (len(wna(x[]))))
setorderv(setDT(TDF),c("rid", "Nna"), order = c(-1,1))
TDF = unique(TDF, incomparables = FALSE, fromLast = FALSE, by="rid")
TDF <- merge(TDF,Rename1Col(IDDT5,"RID","rid")[,c("rid","ND")], by = "rid", all.x = T, all.y = F)

setDF(TDF)
TDF[which(TDF$rid == "rv20_v2_563"),TempRnkcls] = NA
TDF[which(TDF$rid == "rv20_v2_563"),c("phylum","class","order","family")] = c("Lenarviricota","Miaviricetes","Ourlivirales","Botourmiaviridae")

TDF[which(TDF$rid %in% c("rv20_v2_2860","rv20_v2_1282")),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_2860","rv20_v2_1282")),c("phylum","class")] = data.frame(phylum = rep("Pisuviricota",2), class = rep("Pisoniviricetes",2))

TDF[which(TDF$rid %in% c("rv20_v2_214","rv20_v2_565")),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_214","rv20_v2_565")),c("phylum","class")] = data.frame(phylum = rep("Pisuviricota",2), class = rep("Duplopiviricetes",2))

TDF[which(TDF$rid %in% unique(c("ya2_JAAOEH010001907_1","rv20_v2_2491","rv20_v2_3578","rv20_v2_2572","rv20_v2_3579","rv20_v2_6855","rv20_v2_6834","rv20_v2_4230","rv20_v2_5756","rv20_v2_5714","rv20_v2_2491","rv20_v2_2572","Rv4_167678","ya2_JAAOEH010000945_1","Rv4_226693","Rv4_222558"))),TempRnkcls] = NA

TDF[which(TDF$rid %in% c("rv20_v2_7700","rv20_v2_7701","Rv4_188496")),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_7700","rv20_v2_7701","Rv4_188496")),c("phylum","class")] = data.frame(phylum = rep("Pisuviricota",3), class = rep("Stelpaviricetes",3))

TDF[which(TDF$rid %in% c("rv20_v2_4737","rv20_v2_7753","rv20_v2_402")),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_4737","rv20_v2_7753","rv20_v2_402")),c("phylum")] = data.frame(phylum = rep("Pisuviricota",3))

TDF[which(TDF$rid == "rv20_v2_6054"),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_6054")),c("phylum")] = data.frame(phylum = rep("Pisuviricota",1))

TDF[which(TDF$rid %in% c("rv20_v2_1160","rv20_v2_1515","rv20_v2_2064","rv20_v2_6893","rv20_v2_76","rv20_v2_816")),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_1160","rv20_v2_1515","rv20_v2_2064","rv20_v2_6893","rv20_v2_76","rv20_v2_816")),c("phylum")] = c("Negarnaviricota")

TDF[which(TDF$rid %in% c("rv20_v2_3260","rv20_v2_3261","rv20_v2_6505","rv20_v2_3059")),TempRnkcls] = NA
TDF[which(TDF$rid %in% c("rv20_v2_3260","rv20_v2_3261","rv20_v2_6505","rv20_v2_3059")),c("phylum")] = data.frame(phylum = rep("Kitrinoviricota",4))

TDF[which(TDF$rid == "rv20_v2_7187"),TempRnkcls] = NA
TDF[which(TDF$rid == "rv20_v2_7187"),c("phylum","subphylum","class")] = c("Negarnaviricota","Polyploviricotina","Ellioviricetes")


TDF[which(TDF$rid == "ya2_JAAOEH010000746_1"),TempRnkcls] = NA
TDF[which(TDF$rid == "ya2_JAAOEH010000746_1"),c("phylum","class","order","family")] = c("Kitrinoviricota","Tolucaviricetes","Tolivirales","Tombusviridae")

setDT(TDF)
tempTFD = RemoveNACols(data.frame(TDF))
colnames(tempTFD) = BBmisc::capitalizeStrings(colnames(tempTFD))
WriteXlsx(tempTFD,"TDF.xlsx")

tempDCB = RemoveNACols(data.frame(DCB))

colnames(tempDCB) = BBmisc::capitalizeStrings(colnames(tempDCB))
WriteXlsx(tempDCB,"DCB.xlsx")




colnames(DCB) = tolower(colnames(DCB))
colnames(TDF) = tolower(colnames(TDF))

ShCl = SharedCols(DCB,TDF)
DCB2 = plyr::rbind.fill(DCB[,..ShCl], filter(data.frame(TDF), !(nd %in% DCB$nd))[,ShCl])
colnames(DCB2) = BBmisc::capitalizeStrings(colnames(DCB2))
DCB2 <- Rename1Col(DCB2,"Rid","RID")
DCB2 <- Rename1Col(DCB2,"Nd","ND")
# Manually correcting flagged sequences:
AllRnkcls = BBmisc::capitalizeStrings(tolower(RankDF$Rank))
TempRnkcls = intersect(AllRnkcls,colnames(DCB2))

DCB2[which(DCB2$RID == "rv20_v2_563"),TempRnkcls] = NA
DCB2[which(DCB2$RID == "rv20_v2_563"),c("Phylum","Class","Order","Family")] = c("Lenarviricota","Miaviricetes","Ourlivirales","Botourmiaviridae")

DCB2[which(DCB2$RID %in% c("rv20_v2_2860","rv20_v2_1282")),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_2860","rv20_v2_1282")),c("Phylum","Class")] = data.frame(Phylum = rep("Pisuviricota",2), Class = rep("Pisoniviricetes",2))

DCB2[which(DCB2$RID %in% c("rv20_v2_214","rv20_v2_565")),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_214","rv20_v2_565")),c("Phylum","Class")] = data.frame(Phylum = rep("Pisuviricota",2), Class = rep("Duplopiviricetes",2))

DCB2[which(DCB2$RID %in% unique(c("ya2_JAAOEH010001907_1","rv20_v2_2491","rv20_v2_3578","rv20_v2_2572","rv20_v2_3579","rv20_v2_6855","rv20_v2_6834","rv20_v2_4230","rv20_v2_5756","rv20_v2_5714","rv20_v2_2491","rv20_v2_2572","Rv4_167678","ya2_JAAOEH010000945_1","Rv4_226693","Rv4_222558"))),TempRnkcls] = NA

DCB2[which(DCB2$RID %in% c("rv20_v2_7700","rv20_v2_7701","Rv4_188496")),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_7700","rv20_v2_7701","Rv4_188496")),c("Phylum","Class")] = data.frame(Phylum = rep("Pisuviricota",3), Class = rep("Stelpaviricetes",3))

DCB2[which(DCB2$RID %in% c("rv20_v2_4737","rv20_v2_7753","rv20_v2_402")),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_4737","rv20_v2_7753","rv20_v2_402")),c("Phylum")] = data.frame(Phylum = rep("Pisuviricota",3))

DCB2[which(DCB2$RID == "rv20_v2_6054"),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_6054")),c("Phylum")] = data.frame(Phylum = rep("Pisuviricota",1))

DCB2[which(DCB2$RID %in% c("rv20_v2_1160","rv20_v2_1515","rv20_v2_2064","rv20_v2_6893","rv20_v2_76","rv20_v2_816")),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_1160","rv20_v2_1515","rv20_v2_2064","rv20_v2_6893","rv20_v2_76","rv20_v2_816")),c("Phylum")] = c("Negarnaviricota")

DCB2[which(DCB2$RID %in% c("rv20_v2_3260","rv20_v2_3261","rv20_v2_6505","rv20_v2_3059")),TempRnkcls] = NA
DCB2[which(DCB2$RID %in% c("rv20_v2_3260","rv20_v2_3261","rv20_v2_6505","rv20_v2_3059")),c("Phylum")] = data.frame(Phylum = rep("Kitrinoviricota",4))

DCB2[which(DCB2$RID == "rv20_v2_7187"),TempRnkcls] = NA
DCB2[which(DCB2$RID == "rv20_v2_7187"),c("Phylum","Subphylum","Class")] = c("Negarnaviricota","Polyploviricotina","Ellioviricetes")


DCB2[which(DCB2$RID == "ya2_JAAOEH010000746_1"),TempRnkcls] = NA
DCB2[which(DCB2$RID == "ya2_JAAOEH010000746_1"),c("Phylum","Class","Order","Family")] = c("Kitrinoviricota","Tolucaviricetes","Tolivirales","Tombusviridae")


saveRDS(DCB2,"DCB2.RDS")
WriteXlsx(DCB2,"DCB2.xlsx")

{
  TmpWorkTree = WorkTree
  TmpWorkTree = drop.tip(TmpWorkTree,tip = setdiff(TmpWorkTree$tip.label,DCB2$RID))
  ggTmpWorkTree = ggtree(TmpWorkTree, layout = "ape", yscale = "none", ladderize = T)
  info  <- Rename1Col(DCB2, "RID", "label")
  # info <- tmpdf
  pdf("ggTmpWorkTree.19102021.pdf", width = 70, height = 70,title = "DCB pruned w/ RAW taxonomy")
  p2 <- ggTmpWorkTree
  p2$data <- merge(p2$data, info, by = "label", all.x = T, all.y = F)
  p3 = p2 +  geom_tree(aes(color = Order),layout = "ape")
  p3
  dev.off()
}

## VVD's TaxChecks : 
VVDTxCk = ReadXlsx("../Valerian/tree.taxcheckVVD.xlsx")
intersect(VVDTxCk$RID,taxiddf_lin1$RID)

VVDTxCk = merge(VVDTxCk,RemoveNACols(Rename1Col(TDF,"rid","RID")),by="RID",all.x=T,all.y=F)
VVDTxCk <- Rename1Col(VVDTxCk,"nd","ND")

VVDTxCk$CompGbID = (VVDTxCk$GbID == VVDTxCk$gbid)
intersect(DCB$rid, VVDTxCk$RID[wana(VVDTxCk$auto)])


# setdiff(DCB2$ND,IDDT5$ND) #character(0)
IDDT6 <- merge(IDDT5, `dropcols<-`(DCB2,"RID"), by = "ND", all.x = T, all.y = F)
IDDT6 <- Rename1Col(IDDT6,"New_node.label","node.label")
IDDT6$Realm[grep(x = IDDT6$RID, pattern = "al2_", fixed = T)] = "Viruses"
IDDT6$Realm[grep(x = IDDT6$RID, pattern = "ap2_", fixed = T)] = "Viruses"
IDDT6$Phylum[grep(x = IDDT6$RID, pattern = "al2_", fixed = T)] = "Lenarviricota"
IDDT6$Phylum[grep(x = IDDT6$RID, pattern = "ap2_", fixed = T)] = "Pisuviricota"
IDDT6$Class[grep(x = IDDT6$RID, pattern = "ap2_", fixed = T)] = "Stelpaviricetes"
IDDT6$Order[grep(x = IDDT6$RID, pattern = "ap2_", fixed = T)] = "Stellavirales"

saveRDS(IDDT6,"IDDT6.RDS")
WriteXlsx(IDDT6,"IDDT6.xlsx")
fwrite(sep = '\t',file = "IDDT6.tbl",x = IDDT6,quote = T,verbose = T)

ID3xT = merge(IDDT6,data.frame(IDFT)[,c("ND",setdiff(colnames(IDFT),SharedCols(IDFT,IDDT6)))],by = "ND",all.x=T,all.y=F)
fwrite(sep = '\t',file = "ID3xT.tbl",x = ID3xT,quote = T,verbose = T)


# tmpdf = distinct(IDDT6[,c("node",BBmisc::capitalizeStrings(setdiff(ShCl,c("rid","nd"))))])
tmpdf = distinct(IDDT6[,c("node","Phylum","Class","Order","Family","Subfamily")])
tmpdf = tmpdf[wana(tmpdf$node),]

tmptmpdf = plyr::count(tmpdf,"node")
tmpdf$Nna = 0
tmpdf$Nna = apply(tmpdf, MARGIN = 1, FUN = function(x) (length(wna(x[]))))
setorderv(setDT(tmpdf),c("node", "Nna"), order = c(-1,1))
tmpdf = unique(tmpdf, incomparables = FALSE, fromLast = FALSE, by="node")
whd(tmpdf$node)
tmptmpdf = plyr::count(tmpdf,"node")


{
  ggWorkTree = ggtree(WorkTree, layout = "ape", yscale = "none", ladderize = T)
  # info  <- Rename1Col(DCB2, "RID", "label")
  info <- tmpdf
  pdf("ggWorkTree.19102021.pdf", width = 70, height = 70)
  p2 <- ggWorkTree
  p2$data <- merge(p2$data, info, by = "node", all.x = T, all.y = F)
  p3 = p2 +  geom_tree(aes(color = Phylum),layout = "ape")
  p3
  dev.off()
}

# taxonomizr::prepareDatabase(sqlFile = "/media/HDD1/uri/DBs/NCBI_Tax/accessionTaxa.sql",tmpDir = "./tmp/",vocal = T,getAccessions = T)


#### Tree works  ####
{
  WorkTree = as.phylo(read.tree("./210822/ali.210822m/210822m.77510.09.lab_mod.tre"))
  WorkTree$node.label <- p0("node_", seq(1:Nnode(WorkTree, internal.only = T)))
  Nodetbl = data.frame(cbind(WorkTree$edge, WorkTree$edge.length))
  colnames(Nodetbl) = c("Parent_node", "node", "Edge_Length2Parent")
  Nodetbl[, "node.label"] = nodelab(WorkTree, Nodetbl$node);
  Nodetbl[, "IsTip"] = isTip(WorkTree, Nodetbl$node);
  Nodetbl = Nodetbl %>% group_by(Parent_node)
  Nodetbl$DistrFromRoot = get_all_distances_to_root(WorkTree)[Nodetbl$node]
  setorderv(setDT(Nodetbl),c("node"),order = c(-1))
} # Same as in the above (should be versioned as IDDT5+), no need to rerun.

setdiff(Nodetbl$node,IDDT6$node)
Nodetbl = merge(Nodetbl, filter(`dropcols<-`(IDDT6,c("DistrFromRoot","node")),RID %in% WorkTree$tip.label),by = "node.label",all.x=T,all.y=F)



RankDF$Rank = BBmisc::capitalizeStrings(tolower(RankDF$Rank))

Rnkcols = setdiff(intersect(colnames(Nodetbl), (RankDF$Rank)), c("Species", "Subgenus"))
tmpvar = data.frame("Rank" = Rnkcols)
tmpvar$mems = NA
for (i in 1:nrow(tmpvar)) {
  tmpvar$mems[i] = list(unique(data.frame(Nodetbl)[, tmpvar$Rank[i]]))
  print(list(unique(Nodetbl[, tmpvar$Rank[i]])))
}

TaxDT = data.table("Taxa" = na.omit(setdiff(unique(unlist(
  unlist(tmpvar$mems[])
)), "NF")))
TaxDT[, c("LCA", "M_t", "N_t")] = 2.1 # Numeric place holder
TaxDT[, c("LCA_label", "Rank")] = "ASD"  # Char place holder
TaxDT[, c("IsMonophyletic", "IsTip", "IsMonophyletic_wNF")] = F  # Bool place holder
TaxDT[, c("WantedKids")] = vector('list', nrow(TaxDT))   # List place holder 1
TaxDT[, c("UnwantedKids")] = vector('list', nrow(TaxDT)) # List place holder 2
Nodetbl = setDT(Nodetbl)
for (i in 1:nrow(TaxDT)) {
  ttttt = len(tmpvar$Rank[grep(pattern = TaxDT$Taxa[i],
                               x  = tmpvar$mems,
                               fixed = T)])
  if(ttttt != 1){
    print(p0(TaxDT$Taxa[i], " Occures in more than 1 type of rank"))
    next
  }
  
  TaxDT[i, Rank :=  tmpvar$Rank[grep(pattern = TaxDT$Taxa[i],
                                     x  = tmpvar$mems,
                                     fixed = T)]]
  tmprank = which(colnames(Nodetbl) == TaxDT[i, Rank])
  AllanoY = which((Nodetbl[,..tmprank]) == (TaxDT$Taxa[i])) # All node/tip *indices in the table*
  
  if (len(AllanoY) == 0) {
    next
  }
  Allanox = unique(Nodetbl$node.label[AllanoY]) # All node/tip *labels*
  Allano = unique(NLabel2NID(WorkTree, Allanox)) # All node/tip *IDs* (on the tree)
  AllanoTipO = Nodetbl$node.label[AllanoY][which(Nodetbl$IsTip[AllanoY])] # Only tip *IDs*
  if (len(Allano) == 0){
    next
  }
  if (anyNA(Allano)) {
    Allano = Allano[wana(Allano)]
    if (len(Allano) == 0) {
      next
    }
  }
  
  LCA1 = castor::get_mrca_of_set(tree = WorkTree, AllanoTipO)
  
  TaxDT[i, LCA := (LCA1)]
  TaxDT[i, LCA_label :=  NID2NLabel(WorkTree, LCA1)]
  TaxDT[i, IsMonophyletic :=  is_monophyletic(WorkTree, AllanoTipO)]
  TaxDT[i, IsMonophyletic_wNF :=    TaxDT[i, IsMonophyletic]]
  if (  TaxDT[i, IsMonophyletic] == F){
    TaxDT[i, UnwantedKids :=  list(NID2NLabel(WorkTree, setdiff(unique(unlist(
      unlist(GetKids(
        tree =  WorkTree,
        node = LCA1,
        type = "tips"
      ))
    )), Allano)))]
  }
  TaxDT[i, WantedKids :=  list(AllanoTipO)]
  
  tmpppp = try(AnyFalse(unname(unlist(filter(Nodetbl, node.label %in% unlist(TaxDT[i, UnwantedKids]))[, ..tmprank])) %in% c(NA,"NF", TaxDT$Taxa[i])))
  TaxDT[i, IsMonophyletic_wNF :=  !tmpppp]
  TaxDT[i, IsTip :=  isTip(WorkTree, LCA1)]
  TaxDT[i, M_t :=  sum(na.omit(Nodetbl$DistrFromRoot[AllanoY]))] # Phylogenetic distance (cumulative branch length) of the root to each tip *affiliated as this taxa*
  TaxDT[i, N_t :=  length(unlist(TaxDT$WantedKids[i]))]
  if (TaxDT[i, IsMonophyletic_wNF]) {
    # AllkIDS = (NID2NLabel(WorkTree, unique(unlist(
    #   unlist(GetKids(
    #     tree =  WorkTree,
    #     node = LCA1,
    #     type = "tips"
    #   ))
    # ))))
    AllkIDS = union(unlist(TaxDT$WantedKids[i]),unlist(TaxDT$UnwantedKids[i]))
    TaxDT[i, M_t :=  sum(na.omit(filter(Nodetbl, node.label %in% AllkIDS)$DistrFromRoot))] # Phylogenetic distance (cumulative branch length) of the root to each tip *affiliated as this taxa*
    TaxDT[i, N_t :=  len(AllkIDS)]
  }
  print(p0(LCA1, "      -       ", TaxDT$Taxa[i]))
}

TaxDT[which(IsMonophyletic), IsMonophyletic_wNF :=  T]

TaxDT$DistrFromRoot = get_all_distances_to_root(WorkTree)[TaxDT$LCA]
View(TaxDT[, c(
  "Taxa",
  "Rank",
  "LCA",
  "DistrFromRoot",
  "M_t",
  "IsMonophyletic",
  "IsMonophyletic_wNF"
)])

Row2IterOn = TaxDT[, which(IsMonophyletic_wNF)]
# Populate Nodetbl with 2nd strictest LSA definition
for (i in Row2IterOn) {
  tmprank = which(colnames(Nodetbl) == TaxDT[i, Rank])
  tmpkids = unique(unlist(GetKids(WorkTree, TaxDT$LCA[i], type = "all")))
  tmplabs = NID2NLabel(WorkTree, tmpkids)
  Nodetbl[node.label %in% tmplabs, tmprank]  = as.char(TaxDT[i, Taxa])
}
