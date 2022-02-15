# Information  --------------------------------------------------------------------
## Script name: Annotations.r
## Author: Uri Neri
## Email: uri.neri@gmail.com
# Body  --------------------------------------------------------------------
source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
RVMT = "/media/HDD1/uri/RNA_Vir_MTs/RVMT/"
RNA_Vir_MTs_dir = "/media/HDD1/uri/RNA_Vir_MTs/V3/"
versiza = "Annotations"
current_path = p0(RNA_Vir_MTs_dir, versiza, "/")
setwd(current_path)
library(GenomicRanges)
library(plyranges)
library(GenomicTuples)

##### Set Env / Read Info  (local) ####
THREADS = 11
Memory = 11800

IDFT = fread0("../Vfin/Metadata/IDFT.06122021.tsv")
IDFT.fasta = readDNAStringSet("/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/IDFT.fasta")
# IDFT_6frx = readAAStringSet("/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/ORFs/ALL_nuc_3007_6frx.faa")
# names(IDFT_6frx) = HeaderBreaker(data.frame("old_names" = names(IDFT_6frx)),sep  =  " ",clb  =  "old_names")$splitted$X1
AllORFsInfo = fread("../Vfin/ORFs/AllORFsInfo.tsv", header = T , sep = "\t",colClasses=c("partial"="character"))
AllORFsInfo = readRDS("../Vfin/ORFs/ORFsDF.RDS")
AllORFsInfo = AllORFsInfo[which(AllORFsInfo$seqid %in% IDFT$ND),]
AllORFs = readAAStringSet("../Vfin/ORFs/All3007.faa")
AllORFs = AllORFs[intersect(names(AllORFs),AllORFsInfo$ORFID)]
ORFsInfo = distinct(AllORFsInfo[,c("seqid","start","end","strand","ORFID")])
O2N_ORFdf <- distinct(AllORFsInfo[,c("ORFID","OLD_ORFID")])
O2N_ORFdf = distinct(rbind(O2N_ORFdf,fread0("../Vfin/Metadata/Old2New_IDs/ORF.tsv")))
O2N_NUCdf = fread0("../Vfin/Metadata/Old2New_IDs/Contigs.tsv")
# O2N_NC90df = fread0("../Vfin/Metadata/Old2New_IDs/Nucleic_Clusters_90.tsv")
# O2N_NC95df = fread0("../Vfin/Metadata/Old2New_IDs/Nucleic_Clusters_95.tsv")
# O2N_RIDdf = fread0("../Vfin/Metadata/Old2New_IDs/RdRp_IDs.tsv")

hmmsearchcols=c("target_name", "target_accession", "qL", "query_name", "query_accession", "pL", "E-value_fullseq", "score_fullseq", "bias_fullseq", "#", "of", "c-Evalue", "i-Evalue", "bias_thisdomain", "score_thisdomain", "p1", "p2", "q1", "q2", "env_from", "env_to", "acc", "description_of_target")

PolyProts = unique(c("YP_009336655-1243-1333","1ra6A03","4lq3A04","3uqsA04","3h5xA04","5y6zA03","cluster_003204","cluster_004570","YP_009182187-87-805","YP_009330274-238-1037","PF00946","PF00946.21","NP_694468-288-1197","YP_003934933-309-1461","NP_690839-1-875","YP_009094051-171-1093","APG79216-246-1777","P35942-181-1148","AJG39083-75-1055","YP_006390636-200-1050","YP_009336483-1873-2859","YP_002905337-3-1018","YP_009272911-65-1157","ANW72256-210-1321","YP_009305102-50-1080","AIY25916-39-1071","NP_044598-16-1180","NP_066251-10-1089","APG79216246-1777","PR00918","cluster_002482","ALD89111-54-1194","YP_052968-6-1302","YP_004226522-4-1259","AKN56890-1-1262","PF06317.13","APG79225-717-2090","ARF07019-184-2261","AKN56871-186-1869","AJZ68872-213-1628","YP_009362029-223-1914","YP_003620396-1-1220","cluster_003593","NP_694872-178-2213","NP_049362-255-1973","YP_009028573-1209-2746","AKN56888-187-1825","YP_009336824-5-1646","YP_003104764-245-1674","APG79225-717-2090","AGA82737-223-1557","YP_052968-6-1302","ALD89111-54-1194"))
PolyProts = unique(c(PolyProts,scan("Polyprots.txt",what = "\n")))
# 
# CM13 = ReadXlsx("CM13.xlsx",data.table = F)
# CM12 = ReadXlsx("CM12.xlsx",data.table = F)
# CM12 <- rbind.fill(filter(CM12[,SharedCols(CM12,CM13)], !(profile_accession %in% CM13$profile_accession)),CM13)
# CM12$Classified[wana(CM12$New_Name)]  = 1
# CM12$Classified[wna(CM12$New_Name)]  = -1
# 
# CM12$PolyPort_Problematic[CM12$profile_accession %in% PolyProts] <- T
# CM12$PolyPort_Problematic[wna(CM12$PolyPort_Problematic)] = F
# CM12$PolyPort_Problematic <- as.logical(CM12$PolyPort_Problematic)
# CM12$Classified[which(CM12$PolyPort_Problematic)]  = -2
# CM12 <- distinct(CM12)
# ProfCats = data.frame(distinct(CM12[,c("Classified","profile_accession","PolyPort_Problematic","name1")]))
ProfCats = data.frame(distinct(merge(NeoCM3,CM12[,c("name1","profile_accession")],by = "profile_accession",all.x = T,all.y = F)[,c("Classified","profile_accession","PolyPort_Problematic","name1")]))

ProfCats$PolyPort_Problematic = as.logical(ProfCats$PolyPort_Problematic)
ProfCats$Classified[which(ProfCats[,"PolyPort_Problematic"])]  = -2

########## HMMsearch RnaVirdb200205 VS All3007.faa  ########## 
HRVA = (fread("Search_Outputs/awked/awked_domtblout_rnaVirDB_202005_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols,nThread = THREADS))
HRVA = HRVA[wh(HRVA$`E-value_fullseq` <= 1e-4),]
gc()

HRVA = `dropcols<-`(HRVA,c("target_accession", "c-Evalue", "i-Evalue", "score_thisdomain", "env_from", "env_to", "acc", "bias_thisdomain", "description_of_target"))
HRVA = CalcPcoverage(HRVA)
HRVA$ali_Qcov_len = HRVA$q2 - HRVA$q1

colnames(HRVA) = c("OLD_ORFID", "qL", "profile_accession","profile", "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2", "pCoverage","ali_Qcov_len")
HRVA = merge(HRVA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
HRVA = filter(HRVA,ORFID %in% ORFsInfo$ORFID)
HRVA = filter(HRVA,!(profile_accession %in% PolyProts))
# View(distinct(HRVA[,c("profile_accession","pL")]))
HRVA = `dropcols<-`(HRVA,c("r1","r2","r3","profile","OLD_ORFID"))

HRVA = merge(HRVA,ProfCats,by="profile_accession",all.x=T,all.y=F)
which.na(HRVA$Classified)

HRVA = data.table(HRVA)[evalue <= 1e-5,]
HRVA = HRVA[ali_Qcov_len >= 25,]
HRVA = HRVA[score >= 10,]
HRVA = filter(HRVA,!(PolyPort_Problematic))

HRVA = distinct(HRVA)
gc()
framedf = HRVA[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue")] 
colnames(framedf) = c("seqnames","Classified","profile_accession", "start", "end", "qL", "p1", "p2", "pL", "pCoverage", "score", "ali_Qcov_len", "evalue")

framedf$matchID = p0("RVD25.",1:nrow(framedf)) # Not by row names due to data.table. TODO: add input file specific prefix / hash.

setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))
framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}


RDBms = framedf_asgr
RDBms$Analysis_Type = "RNAVirDB200205"
Rna200205df = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(RDBms,stringsAsFactors = F)),factors.as.char  =  T))
rm(RDBms,TMP1,w1,HRVA,framedf_asgr,framedf)
gc()

########## HMMsearch Pfam34a VS All3007.faa  ########## 
HPVA = (fread("Search_Outputs/awked/awked_domtblout_Pfam34A_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
HPVA = HPVA[which(HPVA$`E-value_fullseq` <= 1e-5),]

HPVA = `dropcols<-`(HPVA,c("target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
HPVA = CalcPcoverage(HPVA)
HPVA$ali_Qcov_len = HPVA$q2 - HPVA$q1

colnames(HPVA) = c("OLD_ORFID", "qL", "profile","profile_accession", "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2", "pCoverage","ali_Qcov_len")
HPVA = merge(HPVA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
HPVA = filter(HPVA,ORFID %in% ORFsInfo$ORFID)
HPVA = filter(HPVA,!(profile_accession %in% PolyProts))

NeoCM3 = data.frame(distinct(merge(NeoCM3,HPVA[,c("profile","profile_accession")],by = "profile_accession",all.x = T,all.y = F)))

HPVA = `dropcols<-`(HPVA,c("r1","r2","r3","profile","OLD_ORFID"))

HPVA = merge(HPVA,ProfCats,by="profile_accession", all.x=T, all.y=F, allow.cartesian=T)
w1 = which.na(HPVA$Classified)
if(length(w1)!=0){HPVA  = HPVA[-w1,]}


HPVA = data.table(HPVA)[evalue <= 1e-6,]
HPVA = HPVA[ali_Qcov_len >= 15,]
HPVA = HPVA[score >= 10,]
HPVA = distinct(HPVA)
HPVA = filter(HPVA,!(PolyPort_Problematic))

framedf = HPVA[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue")] 
colnames(framedf) = c("seqnames","Classified","profile_accession", "start", "end", "qL", "p1", "p2", "pL", "pCoverage", "score", "ali_Qcov_len", "evalue")

setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))

framedf$matchID = p0("P34.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.
framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

PfDBms = framedf_asgr
PfDBms$Analysis_Type = "Pfam34A"
PfamDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(PfDBms,stringsAsFactors = F)),factors.as.char  =  T))
# TempDF4Mart = filter(PfamDF,profile_accession == "PF03118.17")
# WriteXlsx(TempDF4Mart,"TempDF4Mart.xlsx")
rm(PfDBms,TMP1,w1,HPVA)
gc()

########## hpc_HMMsearch CDD3.19 VS All3007.faa  ########## 
# SCOPe  >>  ECOD  >> cath >> Pfam34a
CDVA = (fread("Search_Outputs/awked/awked_domtblout_cdd_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
CDVA = CDVA[which(CDVA$`E-value_fullseq` <= 1e-5),]

CDVA = `dropcols<-`(CDVA,c("target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
CDVA = CalcPcoverage(CDVA)
CDVA$ali_Qcov_len = CDVA$q2 - CDVA$q1


colnames(CDVA) = c("OLD_ORFID", "qL", "profile_accession","profile", "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2", "pCoverage","ali_Qcov_len")
CDVA = merge(CDVA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
CDVA = filter(CDVA,ORFID %in% ORFsInfo$ORFID)
CDVA = filter(CDVA,!(profile_accession %in% PolyProts))

CDVA = `dropcols<-`(CDVA,c("r1","r2","r3","profile","OLD_ORFID"))

CDVA = merge(CDVA,ProfCats,by="profile_accession",all.x=T,all.y=F)
CDVA = filter(CDVA,!(PolyPort_Problematic))
w1 = which.na(CDVA$Classified)
if(length(w1)!=0){CDVA  = CDVA[-w1,]}

CDVA = data.table(CDVA)[evalue < 1e-5,]
CDVA = CDVA[ali_Qcov_len >= 15,]
CDVA = CDVA[score >= 12,]
CDVA = distinct(CDVA)

framedf = CDVA[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue")] 

colnames(framedf) = c("seqnames","Classified","profile_accession", "start", "end", "qL", "p1", "p2", "pL", "pCoverage", "score", "ali_Qcov_len", "evalue")

setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))
framedf$matchID = p0("CDD3.19",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.

framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

CDDms = framedf_asgr
CDDms$Analysis_Type = "CDD3.19"
CDDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(CDDms,stringsAsFactors = F)),factors.as.char  =  T))
gc()

########## HMMsearch LysDB VS All3007.faa  ########## 
# LVA = (fread("Search_Outputs/awked/awked_domtblout_LysDB_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
LVA = (fread("All_lys/lys/New_ent2/Merge/awked_domtblout_LysDB_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
LVA = `dropcols<-`(LVA,c("target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
LVA = distinct(LVA)
LVA = CalcPcoverage(LVA)
LVA$ali_Qcov_len = LVA$q2 - LVA$q1

colnames(LVA) = c("ORFID", "qL","profile_accession", "profile",  "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2", "pCoverage","ali_Qcov_len")
LVA$profile = 'Lysis-LysDB'
LVA = `dropcols<-`(LVA,c("r1","r2","r3"))
LVA = distinct(LVA)

LVA = data.table(LVA)[evalue < 0.001,]
LVA = LVA[ali_Qcov_len >= 6,]
LVA = LVA[pCoverage >= 0.01,]

# LVA = merge(LVA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
LVA = filter(LVA,ORFID %in% ORFsInfo$ORFID)

LVA = `dropcols<-`(LVA,c("profile"))
# LVA = merge(LVA,CM12,by="")
LVA = merge(LVA,ProfCats,by="profile_accession",all.x=T,all.y=F)
LVA = filter(LVA,!(PolyPort_Problematic))
which.na(LVA$Classified)

framedf = LVA[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue")] 

colnames(framedf) = c("seqnames", "Classified", "profile_accession","start", "end", "qL","p1", "p2", "pL", "pCoverage", "score", "ali_Qcov_len", "evalue")
framedf$matchID = p0("LysORFs.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.

framedf = framedf[score >= 6,]
setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue","ali_Qcov_len"),order = c(-1,-1,-1,1,-1))
framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F)) 
framedf_asgr = framedf_asgr[unique(TMP1)] 

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = F)) 
framedf_asgr = framedf_asgr[unique(TMP1)] 

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='any', minoverlap = 50, select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}

LysDBms = framedf_asgr[unique(TMP1)] 
LysDBms$Analysis_Type = "Lysis-LysDB"
LysDBms$DB = 8
# LysDBms$Classified = 2

LysDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(LysDBms,stringsAsFactors = F)),factors.as.char  =  T))
gc()

##### MultiSearch - LysDB VS IDFT6Frx.faa #####
psilysout = GenericHitsParser(calc_pcoverage = T,reducecols = T,search_tool = "psiblast",Cull_hits = F,input_was = "ND.frame",breakhdrs = F,Query2Profile = F ,inpt = "/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/All_lys/lys/New_ent2/Out/cated_psiblast_vs_VALL_nuc_3007_sixframe.tsv")
hhlysout = GenericHitsParser(calc_pcoverage = T,reducecols = T,search_tool = "hmmsearch",Cull_hits = F,input_was = "ND.frame",breakhdrs = F,Query2Profile = F ,inpt = "/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/All_lys/lys/New_ent2/Out/awked_domtblout_AllHMMs_hmmsearch_vs_ALL_nuc_3007_sixframe.tsv")
DDlysout = GenericHitsParser(calc_pcoverage = F,reducecols = F,search_tool = "a",colsnms = c("ND.frame",DaimondP_rev),Cull_hits = T,input_was = "ND.frame",breakhdrs = F,Query2Profile = F ,inpt = "/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/All_lys/lys/New_ent2/Out/DiamondP_ALL_nuc_3007_sixframe_vs_singltons.tsv")
DDlysout2 = GenericHitsParser(calc_pcoverage = F,reducecols = F,search_tool = "a",colsnms = c("subject_name","ND.frame","pident","ali_len","p1","p2","q1","q2","pL","qL","mismatch","gapopen","evalue","score"),Cull_hits = T,input_was = "ND.frame",breakhdrs = F,Query2Profile = F ,inpt = "/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/All_lys/lys/New_ent2/Out/DiamondP_singltons_vs_ALL_nuc_3007_sixframe.tsv")
DDllyst = rbind(DDlysout2,DDlysout)

setorderv(setDT(DDllyst),c("ND.frame","score","pident","evalue"),order = c(-1,-1,-1,1))
DDllyst = unique(DDllyst, incomparables=FALSE, fromLast=FALSE,by="ND.frame")
DDllyst = `dropcols<-`(DDllyst,setdiff(colnames(DDllyst),colnames(hhlysout)))
DDllyst = DDllyst[evalue < 0.001 ]
DDllyst$DB=4
LysOut = bind_rows(psilysout,filter(hhlysout,!(ND.frame %in% psilysout$ND.frame)))
LysOut = LysOut[evalue < 0.01 ]
LysOut$DB = 5
LysOut = bind_rows(LysOut,filter(DDllyst,!(ND.frame %in% LysOut$ND.frame)))
LysOut = LysOut[score >= 10  ]

LysOut = HeaderBreakerCb(input_df = data.frame(LysOut),sep = ".",clb ="ND.frame",nclb = c("new_nuc_id","frame") )
LysOut <- merge(LysOut,O2N_NUCdf,by="new_nuc_id")
LysOut$ND.frame <- p0(LysOut$ND,"___",LysOut$frame)
LysOut$frame <- (as.numeric(LysOut$frame))
LysOut$strand = 1
LysOut$strand[wh(LysOut$frame<0)] = -1
LysOut$new_nuc_id <- NULL
LysOut <- merge(LysOut,IDFT[,c("ND","Length")],by="ND",all.x=T,all.y=F)
colnames(LysOut) = c("seqnames","ORFID","profile_accession","qL","pL","p1","p2","Oq1","Oq2","score","evalue","r1","pCoverage","DB","frame","strand","QL")
LysOut <- CalcPcoverage(LysOut)
LysOut <- AAcoor2NAcoor_dFAST(LysOut)
tmpmat <- stringi::stri_split_fixed(str = LysOut$NAcoor, pattern =  "..", simplify = T)
LysOut$NAc2 <- pmax(as.num(tmpmat[,1]),as.num(tmpmat[,2]))
LysOut$NAc1 <- pmin(as.num(tmpmat[,1]),as.num(tmpmat[,2]))

LysOut <- Rename1Col(LysOut,"NAc1","start") 
LysOut <- Rename1Col(LysOut,"NAc2","end")
LysOut$r1 <- NULL

LysOut = merge(LysOut,ProfCats,by="profile_accession",all.x=T,all.y=F)
LysOut = filter(LysOut,!(PolyPort_Problematic))
which.na(LysOut$Classified)

LysOut$strand[wh(LysOut$strand == 1)] = "+"
LysOut$strand[wh(LysOut$strand == -1)] = "-"
framedf = LysOut#[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","score","evalue","DB")] 

framedf$matchID = p0("Lys6frx.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.
setorderv(framedf,cols  =  c("seqnames","Classified","DB","score","evalue"),order = c(-1,-1,-1,-1,1))
framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F)) 
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)] 

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = F)) 
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)] 

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='any', minoverlap = 80, select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
LysDBms2 = framedf_asgr[unique(TMP1)] 
LysDBms2$Analysis_Type = "Lysis-LysDB"
LysDBms2$DB = 5
LysDF2 = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(LysDBms2,stringsAsFactors = F)),factors.as.char  =  T))


LysDF2 = filter(LysDF2,seqnames %in% ORFsInfo$seqid)

# LysDF2$ORFID[which.na(LysDF2$ORFID.x)] = LysDF2$ORFID.y[which.na(LysDF2$ORFID.x)]

LysDF2 = `dropcols<-`(LysDF2,c("ORFID.x","ORFID.y","new_nuc_id","ND"))
###
LysDF3 <- LysDF
colnames(LysDF3) = c("ORFID","q1","q2","width","strand","Classified","profile_accession", "qL", "p1", "p2", "pL", "pCoverage", "score", "ali_Qcov_len", "evalue", "matchID", "Analysis_Type", "DB")
LysDF3 <- `dropcols<-`(LysDF3,c("strand"))


Lys3DF = rbind(LysDF2[,intersect(colnames(LysDF2), colnames(LysDF))], data.frame(LysDF)[,intersect(colnames(LysDF2), colnames(LysDF))])
Lys3DF$DB = 7
# Lys3DF = `droprows<-`(Lys3DF,1121)
# rm(LysDBms,LysDBms2,TMP1,w1,LysDF2,psilysout,DDllyst,DDlysout,DDlysout2,LysOut,hhlysout)
gc()

########## InterProScan All3007.faa VS PHOBIUS, TMHMM and ModdbLite  ########## 
IPR_cols=c("OLD_ORFID", "MD5_Digest", "qL", "Analysis_Type", "profile_accession", "profile", "q1", "q2", "score", "Status", "Date", "IPR_accession", "IPR_annotations")
IPMC = fread("Search_Outputs/IPR/InterProScan_All3007.tsv",sep = "\t",header = F,col.names = IPR_cols)

IPMC = merge(IPMC,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
IPMC = filter(IPMC,ORFID %in% ORFsInfo$ORFID)
IPMC = filter(IPMC,Analysis_Type != "CDD")
IPMC = filter(IPMC,!(profile_accession %in% PolyProts))

IPMC = `dropcols<-`(IPMC,c("qL","Date","MD5_Digest","Status","OLD_ORFID")) #,"i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
IPMC$score[which(IPMC$score =="-")] = 0  
IPMC$evalue = as.numeric(IPMC$score)  
IPMC = `dropcols<-`(IPMC,c("score")) #,"i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))

IPMC$ali_Qcov_len = IPMC$q2 - IPMC$q1

IPMC = data.table(IPMC)[evalue < 1e-6,]

IPMC = distinct(IPMC)

framedf = IPMC

framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ORFID","seqnames")

framedf$matchID = p0("PRINTS",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.
setorderv(framedf,cols  =  c("seqnames","evalue","ali_Qcov_len"),order = c(-1,1,1))
IPMCms = framedf 
IPMCms$DB = 6

IPCCms = IPMCms[which(IPMCms$Analysis_Type %in% c("PRINTS")),]
IPMCms$DB[which(IPMCms$Analysis_Type %in% c("PRINTS"))] = 2.6
IPCCms[,c("score", "p1" ,"p2","pL","width")] = 1
PRINTsdf = setDT(IPCCms)
PRINTsdf = merge(PRINTsdf,ProfCats,by="profile_accession",all.x=T,all.y=F)
PRINTsdf = filter(PRINTsdf,!(PolyPort_Problematic))
which.na(PRINTsdf$Classified)

PRINTsdf$DB = 2.1
PRINTsdf$DB[which(PRINTsdf$Classified == 1 )] = 5

PhobTMModdf = setDT(filter(IPMCms,DB == 6))
PhobTMModdf$matchID = p0("PhobTMMmodd",1:nrow(PhobTMModdf))
rm(IPMC,TMP1,w1,IPCCms,framedf,IPMCms)
gc()


########## HMMsearch CATH VS All3007.faa  ########## 
# After Running Gene3D's cath-resolve-hits and fam_assign.py (inexact names)
HCVA = fread("/media/HDD1/uri/DBs/Gene3d/gene3d_hmmsearch/seqs.crh.csv",sep = ",")
colnames(HCVA) = c(setdiff(colnames(HCVA),"V1"),"r1")
HCVA = `dropcols<-`(HCVA,c("r1","indp-evalue","aligned-regions","boundaries","match-id"))
colnames(HCVA) = c("profile_accession", "cath-superfamily", "OLD_ORFID",  "score","boundaries", "evalue")

HCVA1 = HeaderBreakerCb(input_df = data.frame(HCVA),sep = ",",clb = "boundaries",nclb = c("b1","b2"))

HCVA1 = HeaderBreakerCb(input_df = data.frame(HCVA1),sep = "-",clb = "b1",nclb = c("q1","q4"))
HCVA1 = HeaderBreakerCb(input_df = data.frame(HCVA1),sep = "-",clb = "b2",nclb = c("q3","q2"))

# HCVA1 = HeaderBreakerCb(input_df = data.frame(HCVA),sep = "-",clb = "boundaries",nclb = c("q1","q4","q3","q2"))
HCVA = HCVA1[,c("profile_accession", "cath.superfamily", "OLD_ORFID", "q1","q2", "score", "evalue")]
rm(HCVA1)
HCVA$q1 = as.numeric(HCVA$q1)
HCVA$q2 = as.numeric(HCVA$q2)
HCVA$ali_Qcov_len = HCVA$q2 - HCVA$q1

HCVA = data.table(HCVA)[evalue < 1e-6,]
HCVA = HCVA[ali_Qcov_len >= 10,]
HCVA = HCVA[score >= 8,]
HCVA = merge(HCVA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
HCVA = filter(HCVA,ORFID %in% ORFsInfo$ORFID)
HCVA = filter(HCVA,!(profile_accession %in% PolyProts))

HCVA = `dropcols<-`(HCVA,c("OLD_ORFID"))

HCVA = merge(HCVA,ProfCats,by="profile_accession",all.x=T,all.y=F)
HCVA = filter(HCVA,!(PolyPort_Problematic))
which.na(HCVA$Classified)


framedf = HCVA
colnames(framedf) = c("seqnames","profile_accession", "cath.superfamily","start", "end",  "score","evalue", "ali_Qcov_len","Classified", "ND", "strand")
framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ORFID","seqnames")

Cathdf = framedf
Cathdf$Analysis_Type = "Cath"
Cathdf$DB = 2.4
Cathdf$matchID = p0("CATH.",1:nrow(Cathdf))
Cathdf$DB[which(Cathdf$Classified == 1 )] = 4
colnames(Cathdf) = c("seqnames","profile_accession", "cath.superfamily","start","end","score","evalue","ali_Qcov_len","Classified","ND","strand","Analysis_Type","DB","matchID"  )

rm(CathDBms,TMP1,w1,framedf_asgr,framedf,HCVA)
gc()


########## HMMsearch ecod VS All3007.faa  ########## 
EDVA = (fread("Search_Outputs/awked/awked_domtblout_ecod_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
EDVA = EDVA[which(EDVA$`E-value_fullseq` <= 1e-5),]

EDVA = merge(EDVA,Rename1Col(O2N_ORFdf,"OLD_ORFID","target_name"),by="target_name",all.x=T,all.y=F)

EDVA = `dropcols<-`(EDVA,c("OLD_ORFID","target_name","target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
EDVA = CalcPcoverage(EDVA)
EDVA$ali_Qcov_len = EDVA$q2 - EDVA$q1
colnames(EDVA) = c("qL","profile", "profile_accession", "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2", "ORFID", "pCoverage","ali_Qcov_len")
EDVA = filter(EDVA,ORFID %in% ORFsInfo$ORFID)
EDVA = filter(EDVA,!(profile_accession %in% PolyProts))

EDVA = `dropcols<-`(EDVA,c("r1","r2","r3"))

EDVA = data.table(EDVA)[evalue < 1e-7,]
EDVA = EDVA[ali_Qcov_len >= 10,]
EDVA = EDVA[score >=9,]

EDVA = distinct(EDVA)
EDVA = merge(EDVA,ProfCats,by="profile_accession",all.x=T,all.y=F,allow.cartesian=TRUE)
EDVA = filter(EDVA,!(PolyPort_Problematic))

which.na(EDVA$Classified)

EDVA$DB = 4 
EDVA$DB[which(EDVA$Classified %in% c(1,2,3,4) )] = 5

framedf = EDVA[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue","DB")] 

framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ORFID","seqnames")

framedf$matchID = p0("ECOD.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.
setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))

framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}
EcodDBms = framedf_asgr
EcodDBms$Analysis_Type = "ECOD"
EcodDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(EcodDBms,stringsAsFactors = F)),factors.as.char  =  T))

rm(EcodDBms,TMP1,w1,framedf_asgr,framedf,EDVA)
gc()

########## HMMsearch InterDomainProfiles VS All3007.faa  ########## 
IDPA = (fread("Search_Outputs/awked/awked_domtblout_interdomain_hmmsearch_vs_All3007.tsv", sep = "\t",header = F,col.names = hmmsearchcols))
IDPA = IDPA[which(IDPA$`E-value_fullseq` <= 1e-6),]

IDPA = `dropcols<-`(IDPA,c("target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
IDPA$ali_Qcov_len = IDPA$q2 - IDPA$q1
IDPA = na.omit(IDPA)
IDPA = CalcPcoverage(IDPA)

colnames(IDPA) = c("OLD_ORFID", "qL", "profile_accession", "profile", "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2","ali_Qcov_len", "pCoverage")
IDPA = merge(IDPA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
IDPA = filter(IDPA,ORFID %in% ORFsInfo$ORFID)

IDPA = filter(IDPA,!(profile_accession %in% PolyProts))

IDPA = `dropcols<-`(IDPA,c("OLD_ORFID","profile"))

IDPA = data.table(IDPA)[evalue < 1e-9,]
IDPA = IDPA[ali_Qcov_len >= 15,]
IDPA = IDPA[score >= 10,]

IDPA = IDPA[pCoverage >= 0.05,]
IDPA$profile = "Unclassified-InterDomain"
IDPA = distinct(IDPA)
IDPA = merge(IDPA,ProfCats,by="profile_accession",all.x=T,all.y=F)
IDPA = filter(IDPA,!(PolyPort_Problematic))
which.na(IDPA$Classified)

IDPA$DB = 1
IDPA$DB[which(IDPA$Classified %in% c(1,2,3))] = 5

framedf = IDPA[,c("ORFID","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","DB","evalue")] 

framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ORFID","seqnames")

framedf$matchID = p0("IntDom.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.

setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))

framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

InterDoms = framedf_asgr
InterDoms$Analysis_Type = "InterDomains"
InterDomsdf = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(InterDoms,stringsAsFactors = F)),factors.as.char  =  T))
rm(InterDoms,TMP1,w1,framedf_asgr,framedf,IDPA)
gc()

########## HMMsearch SCOPe(1.75) VS All3007.faa  ########## 
SCVA = (fread("Search_Outputs/awked/awked_domtblout_SCOPe_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
SCVA = SCVA[which(SCVA$`E-value_fullseq` <= 1e-5),]

SCVA = `dropcols<-`(SCVA,c("target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
SCVA = CalcPcoverage(SCVA)
SCVA$ali_Qcov_len = SCVA$q2 - SCVA$q1

colnames(SCVA) = c("ASDASD", "qL","profile", "profile_accession", "pL","evalue", "score", "r1", "r2","r3","p1","p2", "q1", "q2", "pCoverage","ali_Qcov_len")
SCVA$OLD_ORFID = gsub(pattern = ".", replacement = "|", x = SCVA$ASDASD,fixed = T)

SCVA = merge(SCVA,O2N_ORFdf,by="OLD_ORFID",all.x=T,all.y=F)
SCVA = filter(SCVA,ORFID %in% ORFsInfo$ORFID)
SCVA = filter(SCVA,!(profile_accession %in% PolyProts))
SCVA = merge(SCVA,ProfCats,by="profile_accession",all.x=T,all.y=F)
SCVA = filter(SCVA,!(PolyPort_Problematic))
SCVA = `dropcols<-`(SCVA,c("r1","r2","r3","ASDASD","OLD_ORFID","profile","PDBID","Class","Function_ID"))
w1 = which.na(SCVA$Classified)
if(length(w1)!=0){SCVA  = SCVA[-w1,]}

SCVA$DB = 4
SCVA$DB[which(SCVA$Classified %in% c(1,2,3))] = 5

SCVA = data.table(SCVA)[evalue < 1e-6,]
SCVA = SCVA[ali_Qcov_len >= 9,]
SCVA = SCVA[score >= 9,]

SCVA = distinct(SCVA)

framedf = SCVA[,c("ORFID","DB","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue")] 

framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ORFID","seqnames")

framedf$matchID = p0("SCOPe.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.

setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))

framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

SCOPeDBms = framedf_asgr
SCOPeDBms$Analysis_Type = "SCOPe"
SCOPeDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(SCOPeDBms,stringsAsFactors = F)),factors.as.char  =  T))

rm(SCOPeDBms,TMP1,w1,framedf_asgr,framedf,SCVA)
gc()


########## HMMsearch Splitted_Profs VS All3007.faa  ##########  
SPVA = (fread("Search_Outputs/awked/awked_domtblout_SP_hmmsearch_vs_All3007.tsv",sep = "\t",header = F,col.names = hmmsearchcols))
SPVA = SPVA[which(SPVA$`E-value_fullseq` <= 1e-5),]

SPVA = `dropcols<-`(SPVA,c("target_accession","c-Evalue","i-Evalue","score_thisdomain","env_from",	"env_to",	"acc","bias_thisdomain","description_of_target"))
SPVA = CalcPcoverage(SPVA)
SPVA$ali_Qcov_len = SPVA$q2 - SPVA$q1

SPVA$ORFID = SPVA$target_name
SPVA$target_name = NULL
SPVA$profile_accession = SPVA$query_name
SPVA = filter(SPVA,ORFID %in% ORFsInfo$ORFID)
SPVA = filter(SPVA,!(profile_accession %in% PolyProts))

SPVA = `dropcols<-`(SPVA,c("query_name","bias_fullseq","#","query_accession","OLD_ORFID","profile","PDBID","Class","of"))
colnames(SPVA) = c("qL", "pL", "evalue","score","p1", "p2", "q1", "q2","pCoverage", "ali_Qcov_len","ORFID","profile_accession" )
SPVA = merge(SPVA,ProfCats,by="profile_accession",all.x=T,all.y=F)
w1 = which.na(SPVA$Classified)
if(length(w1)!=0){SPVA  = SPVA[-w1,]}

SPVA$DB = 7

SPVA$DB[which(SPVA$Classified %in% c(1,2,3) )] = 6

SPVA = data.table(SPVA)[evalue < 1e-5,]
SPVA = SPVA[ali_Qcov_len >= 6,]
SPVA = SPVA[score >= 9,]

SPVA = distinct(SPVA)


framedf = SPVA[,c("ORFID","DB","Classified","profile_accession","q1","q2","qL","p1","p2","pL","pCoverage","score","ali_Qcov_len","evalue")] 

framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ORFID","seqnames")
framedf$matchID = p0("Splitted_Prof.",1:nrow(framedf)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.
setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))

framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}


SPDBms = framedf_asgr
SPDBms$Analysis_Type = "Splitted_Prof"
SPDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(SPDBms,stringsAsFactors = F)),factors.as.char  =  T))
SPDF$Classified = 1
rm(SPDBms,TMP1,w1,framedf_asgr,framedf,SPVA)
gc()

##### 05-.12.2021 Insert NewReps (+segmented) into AllDF #####
tmpcols <- c("ND.frame","qL","profile_accession","profile_name","pL","evalue","score","bias_fullseq","p1","p2","q1","q2")
NrpsVSNVPCp <- fread(na.strings = "",header = F,fill=T,"awked_domtblout_Mprof_hmmsearch_vs_IDFT6frx.tsv",nThread = THREADS)
colnames(NrpsVSNVPCp) <- tmpcols
# NrpsVSNVPCp <- fread("NVPC/SingleTypeProfs/test_hpchmmsearch/domtblout_Mprof_hmmsearch_vs_NewReps.tsv",na.strings = "",header = F,drop = c(2,25))
# colnames(NrpsVSNVPCp) <- c("ND.frame","qL","profile_accession","profile_name","pL","evalue","score","bias_fullseq","#","of","c-Evalue","i-Evalue","bias_thisdomain","score_thisdomain","p1","p2","q1","q2","env_from","env_to","acc","Desc","R")
# NrpsVSNVPCp <- `dropcols<-`(NrpsVSNVPCp,c("R","env_from","env_to","acc","bias_fullseq","#","of","c-Evalue","i-Evalue","bias_thisdomain"))
NrpsVSNVPCp$score = as.double(NrpsVSNVPCp$score)
NrpsVSNVPCp$p1 = as.num(NrpsVSNVPCp$p1)
NrpsVSNVPCp$p2 = as.num(NrpsVSNVPCp$p2)
NrpsVSNVPCp$q1 = as.num(NrpsVSNVPCp$q1)
NrpsVSNVPCp$q2 = as.num(NrpsVSNVPCp$q2)
NrpsVSNVPCp$evalue = as.double(NrpsVSNVPCp$evalue)
NrpsVSNVPCp$pL = as.num(NrpsVSNVPCp$pL)
NrpsVSNVPCp$qL = as.num(NrpsVSNVPCp$qL)
NrpsVSNVPCp = NrpsVSNVPCp[wana(NrpsVSNVPCp$p1),]
NrpsVSNVPCp = NrpsVSNVPCp[wh(NrpsVSNVPCp$q1 !="--incdomE"),]   

NrpsVSNVPCp$profile_accession[wh(NrpsVSNVPCp$profile_name != "-")] <- NrpsVSNVPCp$profile_name[wh(NrpsVSNVPCp$profile_name != "-")]
NrpsVSNVPCp$profile_name <- NULL

NrpsVSNVPCp$profile_accession <- gsub(pattern = ".faa.msa.Cons.msa",replacement = "",x = NrpsVSNVPCp$profile_accession)
NrpsVSNVPCp$profile_accession <- gsub(pattern = ".msa",replacement = "",x = NrpsVSNVPCp$profile_accession)
PolyProts <- c("cluster_003204", "cluster_004570", "YP_009182187-87-805", "YP_009330274-238-1037", "PF00946", "PF00946.21", "NP_694468-288-1197", "YP_003934933-309-1461", "NP_690839-1-875", "YP_009094051-171-1093", "APG79216-246-1777", "P35942-181-1148", "AJG39083-75-1055", "YP_006390636-200-1050", "YP_009336483-1873-2859", "YP_002905337-3-1018", "YP_009272911-65-1157", "ANW72256-210-1321", "YP_009305102-50-1080", "AIY25916-39-1071", "NP_044598-16-1180", "NP_066251-10-1089", "APG79216246-1777", "PR00918", "cluster_002482", "ALD89111-54-1194", "YP_052968-6-1302", "YP_004226522-4-1259", "AKN56890-1-1262", "PF06317.13", "APG79225-717-2090", "ARF07019-184-2261", "AKN56871-186-1869", "AJZ68872-213-1628", "YP_009362029-223-1914", "YP_003620396-1-1220", "cluster_003593", "NP_694872-178-2213", "NP_049362-255-1973", "YP_009028573-1209-2746", "AKN56888-187-1825", "YP_009336824-5-1646", "YP_003104764-245-1674", "AGA82737-223-1557", "putative_polyprotein", "NP_573541-1-676", "AFN73048-160-772", "AII01805-1323-1950", "YP_004598981-3600-4280", "ARF07019-184", "APG79225", "YP_009330274", "cd21530", "YP_002308505-4545-4591", "cluster_004196", "cluster_002083", "cluster_003439", "cluster_004301", "cluster_002164", "cluster_000913", "cluster_000945", "cluster_002691", "cd21593", "cd21588")
NrpsVSNVPCp = filter(NrpsVSNVPCp,!(profile_accession %in% PolyProts))
NrpsVSNVPCp$ali_Qcov_len = NrpsVSNVPCp$q2 - NrpsVSNVPCp$q1
NrpsVSNVPCp = CalcPcoverage(NrpsVSNVPCp)
NrpsVSNVPCp <- NrpsVSNVPCp[evalue < 0.00001]
NrpsVSNVPCp <- NrpsVSNVPCp[score > 10]
NrpsVSNVPCp <- NrpsVSNVPCp[ali_Qcov_len > 10]
gc()

framedf = NrpsVSNVPCp

framedf <- Rename1Col(framedf,"q1","start")
framedf <- Rename1Col(framedf,"q2","end")
framedf <- Rename1Col(framedf,"ND.frame","seqnames")
setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))

framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = wna(TMP1)
if(len(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()


testtest=as.char(unique(framedf_asgr@seqnames))
set1 <- testtest[1:60000]
set2 <- testtest[60000:120000]
set3 <- testtest[120000:180000]
set4 <- testtest[180000:240000]
set5 <- testtest[240000:300000]
set6 <- testtest[300000:360000]
set7 <- testtest[360000:420000]
set8 <- testtest[420000:480000]
set9 <- testtest[480000:526405]

MastersetLst <- framedf_asgr[0]
for (setn in list(set1,set2,set3,set4,set5,set6,set7,set8,set9)){
  workset=filter(framedf_asgr,seqnames %in% setn)
  TMP1 = (GenomicRanges::findOverlaps(workset,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = wna(TMP1)
  if(len(w1)!=0){TMP1[w1] = w1}
  workset = workset[unique(TMP1)]
  gc()
  for(i in c(1:4)){
    TMP1 = (GenomicRanges::findOverlaps(workset,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
    w1 = wna(TMP1)
    if(len(w1)!=0){TMP1[w1] = w1}
    workset = workset[unique(TMP1)]
    gc()
  }
  MastersetLst <- c(MastersetLst,workset)
  rm(workset,TMP1,w1)
  gc()
}
saveRDS(MastersetLst,"MastersetLst.RDS")


framedf_asgr <- readRDS("MastersetLst.RDS")
for(i in c(1:1)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 50,minoverlap = 50, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = wna(TMP1)
  if(len(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}


for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 50, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = wna(TMP1)
  if(len(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}


framedf_asgr$Analysis_Type = "NVPCpl"
NrpsDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(framedf_asgr,stringsAsFactors = F)),factors.as.char  =  T))

# {NeoCM <- distinct(AllDF2[,SharedCols(AllDF2,CM13)])
# NeoCM <- rbind.fill(data.frame(NeoCM),data.frame(workDF)[wh(!(workDF$profile_accession %in% NeoCM$profile_accession)),c(SharedCols(NVPC_Profs,NeoCM))])
# 
# NeoCM <- rbind.fill(data.frame(NeoCM),data.frame(NVPC_Profs)[wh(!(NVPC_Profs$profile_accession %in% NeoCM$profile_accession)),c(SharedCols(NVPC_Profs,NeoCM))])
# NeoCM <- rbind.fill(data.frame(NeoCM),data.frame(CM14)[wh(!(CM14$profile_accession %in% NeoCM$profile_accession)),c(SharedCols(NVPC_Profs,NeoCM))])
# NeoCM <- rbind.fill(data.frame(NeoCM),data.frame(CM13)[wh(!(CM13$profile_accession %in% NeoCM$profile_accession)),c(SharedCols(NVPC_Profs,NeoCM))])
# NeoCM <- rbind.fill(data.frame(NeoCM),data.frame(CM12)[wh(!(CM12$profile_accession %in% NeoCM$profile_accession)),c(SharedCols(NVPC_Profs,NeoCM))])
# NeoCM <- distinct(NeoCM[,c("New_Name","Comment","profile_accession","PolyPort_Problematic","Classified","ClanID")])
# View(MiniDarProf[wh(MiniDarProf$profile_accession %in% NeoCM[NeoCM$New_Name == NeoCM$profile_accession,"profile_accession"]),])
# NeoCM$New_Name[wh(NeoCM$profile_accession %in% c("AAK69629-611-1065","ABX79997-246-759","YP_009344970-4139-4220"))] <- "RdRp"
# NeoCM$New_Name[wh(NeoCM$profile_accession %in% c("ACB45490-1730-1923"))] <- "HCV_NS4b"
# NeoCM$New_Name[wh(NeoCM$ClanID %in% c("Clan_0604"))] <- "HCV_env"
# 
# NeoCM$New_Name[wh(NeoCM$ClanID %in% c("Clan_1244"))] <- "HCV_NS5A"
# NeoCM$New_Name[wh(NeoCM$ClanID %in% c("Clan_1441"))] <- "VSR_dsRBD"
# NeoCM$New_Name[wh(NeoCM$ClanID %in% c("Clan_1140"))] <- "CP_SJR"
# 
# w1 <- intersect(intersect(wh(NeoCM$New_Name != NeoCM$profile_accession),wana(NeoCM$New_Name)),wh(NeoCM$Classified %in% c(NA,-1)))
# NeoCM$Classified[w1] <- "1"
# }
# NeoCM$Classified[wh(NeoCM$profile_accession %in% c("ACB45490-193-382", "T08839-26-198","AAK69629-611-1065","ABX79997-246-759","ACB45490-1730-1923"))] <- 1


NeoCM3 <- fread("../../RVMT/NeoCM3.tsv",header = T,sep="\t",na.strings = "")

NrpsDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(framedf_asgr,stringsAsFactors = F)),factors.as.char  =  T))
NrpsDF$New_Name=NULL
NrpsDF$Classified=NULL

NrpsDF <- merge(NrpsDF,NeoCM3,by="profile_accession",all.x=T,all.y=F)
NrpsDF$New_Name[grep(fixed=T,"RdRp_set2_cluster_",x =NrpsDF$profile_accession )] = "RdRp"
NrpsDF$New_Name[grep(fixed=T,"set.22",x =NrpsDF$profile_accession )] = "RdRp"
NrpsDF$New_Name[grep(fixed=T,"set.25",x =NrpsDF$profile_accession )] = "RdRp"
NrpsDF$New_Name[grep(fixed=T,"yangshan",x =NrpsDF$profile_accession )] = "RdRp"
NrpsDF$New_Name[wh(NrpsDF$profile_accession =="Lysin_M15")] = "Lysin_M15"
NrpsDF$New_Name[wh(NrpsDF$profile_accession =="Lysin_M34")] = "Lysin_M34"

NrpsDF$Classified[grep(fixed=T,"set.22",x =NrpsDF$profile_accession )] = 3
NrpsDF$Classified[grep(fixed=T,"set.25",x =NrpsDF$profile_accession )] = 3
NrpsDF$Classified[grep(fixed=T,"yangshan",x =NrpsDF$profile_accession )] = 3
NrpsDF$Classified[grep(fixed=T,"RdRp_set2_cluster_",x =NrpsDF$profile_accession )] = 0.1

NrpsDF$Classified[wh(NrpsDF$profile_accession =="Lysin_M15")] = 3
NrpsDF$Classified[wh(NrpsDF$profile_accession =="Lysin_M34")] = 3

# 
# 
# NeoCM2 <- rbind.fill(data.frame(NeoCM),distinct(data.frame(NrpsDF)[wh(!(NrpsDF$profile_accession %in% NeoCM$profile_accession)),c(SharedCols(NrpsDF,NeoCM))]))
# NeoCM3 <- distinct(NeoCM2[,c("New_Name","Comment","profile_accession","PolyPort_Problematic","Classified","ClanID")])
# NeoCM3$New_Name[wh(NeoCM3$profile_accession %in% LysOut$profile_accession)[wna(NeoCM3$New_Name[wh(NeoCM3$profile_accession %in% LysOut$profile_accession)])]] = "Lysis"
# 
# View(NeoCM3[wna(NeoCM3$New_Name),])
# w1 <- intersect(intersect(wh(NrpsDF$New_Name != NrpsDF$profile_accession),wana(NrpsDF$New_Name)),wh(NrpsDF$Classified %in% c(NA,-1)))
# NrpsDF$Classified[w1] <- "1"


# 
# w1 <- intersect(intersect(wh(NeoCM3$New_Name != NeoCM3$profile_accession),wana(NeoCM3$New_Name)),wh(NeoCM3$Classified %in% c(NA,-1)))
# NeoCM3$Classified[w1] <- 1
# NeoCM3$Classified[wh(NeoCM3$Classified == "-1")] <- -1
# NeoCM3$Classified[wna(NeoCM3$New_Name)] <- -2
# # NeoCM3$Classified <- as.num(NeoCM3$Classified)
View(merge(merge(data.frame(plyr::count(NrpsDF[wh(NrpsDF$Classified %in% c(-1,"-1",NA)),c("profile_accession","New_Name")])),DarProf_annotations,by="profile_accession",all.x=T,all.y=F),NeoCM3,by="profile_accession",all.x=T,all.y=F))

WriteXlsx(filepath = "../../DariusPreds_x_Prof_Counts.xslx",dfdt = merge(merge(data.frame(plyr::count(NrpsDF[wh(NrpsDF$Classified %in% c(-1,"-1",NA)),c("profile_accession","New_Name")])),DarProf_annotations,by="profile_accession",all.x=T,all.y=F),NeoCM3,by="profile_accession",all.x=T,all.y=F))


# 
# len(wh(NrpsDF$Classified < 0))

gc()

NrpsDF <- Rename1Col(NrpsDF,"start","q1")
NrpsDF <- Rename1Col(NrpsDF,"end","q2")
NrpsDF <- Rename1Col(NrpsDF,"seqnames","ND.frame")

NrpsDF$ND <- str_split_fixed(string = NrpsDF$ND.frame,pattern = fixed("."),n=2)[,1]
NrpsDF$frame <-   as.num(str_split_fixed(string = NrpsDF$ND.frame,pattern = fixed("."),n=2)[,2])
NrpsDF$strand <- NULL
NrpsDF <- merge(NrpsDF,IDFT[,c("ND","Length")],by="ND",all.x=T,all.y=F)
fDF <- wh(NrpsDF$frame>0)
rDF <- wh(NrpsDF$frame<0)
NrpsDF$strand <- "+"
NrpsDF$strand[rDF] <- "-"
NrpsDF$n1 <- -1
NrpsDF$n2 <- -1

NrpsDF$n1[fDF] = (NrpsDF$frame[fDF] + (NrpsDF$q1[fDF]*3)) - 3
NrpsDF$n2[fDF] = (NrpsDF$frame[fDF] + (NrpsDF$q2[fDF]*3)) - 1

NrpsDF$n1[rDF] = (NrpsDF$Length[rDF] + NrpsDF$frame[rDF] - (NrpsDF$q1[rDF]*3)) + 4
NrpsDF$n2[rDF] = (NrpsDF$Length[rDF] + NrpsDF$frame[rDF] - (NrpsDF$q2[rDF]*3)) + 2
NrpsDF$ORFID <- p0(NrpsDF$ND,"__",NrpsDF$frame,"__",NrpsDF$q1,"..",NrpsDF$q2)
fwrite(NrpsDF,"AllDF.07122021.tsv",sep = "\t")



NrpsDF$n3 <- pmax(NrpsDF$n1,NrpsDF$n2)
NrpsDF$n1 <- pmin(NrpsDF$n1,NrpsDF$n2)
NrpsDF$n2 <- NrpsDF$n3
NrpsDF$n3 <- NULL
NrpsDF <- Rename1Col(NrpsDF,"n1","start") 
NrpsDF <- Rename1Col(NrpsDF,"n2","end")
saveRDS(NrpsDF,"NrpsDF.RDS")


NrpsDF[,c("N1","N2")] = -1
NrpsDF$N1 <- NrpsDF$start
NrpsDF$N2 <- NrpsDF$end
rows2iter <- wh(NrpsDF$strand=="-")
setDT(NrpsDF)
NrpsDF$SelfDup <- F
patternS <-
  reverseComplement(narrow(IDFT.fasta[NrpsDF$ND[rows2iter]],start = unlist(pmin(NrpsDF[rows2iter,c("end")],NrpsDF[rows2iter,c("start")])),end =  unlist(pmax(NrpsDF[rows2iter,c("end")],NrpsDF[rows2iter,c("start")]))))
subjectS <- 
  reverseComplement(IDFT.fasta[NrpsDF$ND[rows2iter]])
TMPestDF <- data.frame(patternS = as.char(patternS), subjectS = as.char(subjectS))
TMPestDF$coor <- str_locate_all(pattern = TMPestDF$patternS,string = TMPestDF$subjectS)
for(ix in 1:len(rows2iter)){
  # rawseq <- IDFT.fasta[NrpsDF$ND[ix]]
  # revcomphit <- reverseComplement(narrow(rawseq,start = min(NrpsDF[ix,c("start","end")]),end = max(NrpsDF[ix,c("start","end")])))
  # ReCorr <- vmatchPattern(pattern = BString(x = as.char(revcomphit)),subject = reverseComplement(rawseq), algorithm="naive-exact")
  if(nrow(TMPestDF$coor[[ix]]) > 1){
    NrpsDF[rows2iter[ix], SelfDup:= T]
    next
  }
  NrpsDF[rows2iter[ix], N1:=  TMPestDF$coor[[ix]][1]]
  NrpsDF[rows2iter[ix], N2:= TMPestDF$coor[[ix]][2]]
  print(rows2iter[ix])
}



#### Add manual annotations ####
p0002MA <- ReadXlsx("ManualHitsFromMart.xlsx")
p0002MA <- merge(p0002MA,ORFsInfo,by="ORFID",all.x=F,all.y=F)
p0002MA$q1 <- 1
p0002MA$q2 <- ((p0002MA$end+1 - p0002MA$start) / 3) - 2
p0002MA$score=NULL
p0002MA$Classified <- 2
p0002MA$PolyPort_Problematic <- F
setDF(p0002MA)
p0002MA <- `dropcols<-`(p0002MA,c("seqid","end","start"))

p0002MA <- Rename1Col(p0002MA,"q1","start")
p0002MA <- Rename1Col(p0002MA,"q2","end")
p0002MA <- Rename1Col(p0002MA,"ORFID","seqnames")
p0002MA$matchID = p0("ManualMart.",1:nrow(p0002MA)) # Not by rownames due to data.table. TODO: add input file specific prefix / hash.
CM12 <- rbind.fill(CM12,distinct(p0002MA[,SharedCols(p0002MA,CM12)]))

##### 22-.11.2021 Insert RdRps into AllDF #####
someDF1 <- fread("../Wolf/someDF1.tbl",na.strings = "")
fDF <- wh(someDF1$frame>0)
rDF <- wh(someDF1$frame<0)
someDF1$NAc1
someDF1$n1 <- -1
someDF1$n2 <- -1

someDF1$n1[fDF] = (someDF1$frame[fDF] + (someDF1$q1[fDF]*3)) - 3
someDF1$n2[fDF] = (someDF1$frame[fDF] + (someDF1$q2[fDF]*3)) - 1

someDF1$n1[rDF] = (someDF1$Length[rDF] + someDF1$frame[rDF] - (someDF1$q1[rDF]*3)) + 4
someDF1$n2[rDF] = (someDF1$Length[rDF] + someDF1$frame[rDF] - (someDF1$q2[rDF]*3)) + 2
someDF1$ORFID <- p0(someDF1$ND,"__",someDF1$frame,"__",someDF1$q1,"..",someDF1$q2)

RdRpDF <- distinct(someDF1[,c("ND","ORFID","n1","n2","profile","evalue","score","p1","p2","pL","strand","frame")])  
RdRpDF=AAcoor2NAcoor_df(dfdt = someDF1,frmcol = "frame",AAc1col = "q1",AAc2col = "q2",outdf = T)
RdRpDF$n3 <- pmax(RdRpDF$n1,RdRpDF$n2)
RdRpDF$n1 <- pmin(RdRpDF$n1,RdRpDF$n2)
RdRpDF$n2 <- RdRpDF$n3
RdRpDF$n3 <- NULL
RdRpDF <- Rename1Col(RdRpDF,"n1","start") 
RdRpDF <- Rename1Col(RdRpDF,"n2","end")
saveRDS(RdRpDF,"RdRpDF.RDS")


RdRpDF[,c("N1","N2")] = -1
RdRpDF$N1 <- RdRpDF$start
RdRpDF$N2 <- RdRpDF$end
rows2iter <- wh(RdRpDF$strand=="-")
setDT(RdRpDF)
RdRpDF$SelfDup <- F
patternS <-
  reverseComplement(narrow(IDFT.fasta[RdRpDF$ND[rows2iter]],start = unlist(pmin(RdRpDF[rows2iter,c("end")],RdRpDF[rows2iter,c("start")])),end =  unlist(pmax(RdRpDF[rows2iter,c("end")],RdRpDF[rows2iter,c("start")]))))
subjectS <- 
  reverseComplement(IDFT.fasta[RdRpDF$ND[rows2iter]])
TMPestDF <- data.frame(patternS = as.char(patternS), subjectS = as.char(subjectS))
TMPestDF$coor <- str_locate_all(pattern = TMPestDF$patternS,string = TMPestDF$subjectS)
for(ix in 1:len(rows2iter)){
  # rawseq <- IDFT.fasta[RdRpDF$ND[ix]]
  # revcomphit <- reverseComplement(narrow(rawseq,start = min(RdRpDF[ix,c("start","end")]),end = max(RdRpDF[ix,c("start","end")])))
  # ReCorr <- vmatchPattern(pattern = BString(x = as.char(revcomphit)),subject = reverseComplement(rawseq), algorithm="naive-exact")
  if(nrow(TMPestDF$coor[[ix]]) > 1){
    RdRpDF[rows2iter[ix], SelfDup:= T]
    next
  }
  RdRpDF[rows2iter[ix], N1:=  TMPestDF$coor[[ix]][1]]
  RdRpDF[rows2iter[ix], N2:= TMPestDF$coor[[ix]][2]]
  print(rows2iter[ix])
}


##### 28.11.2021 Add HHsearch Reps6Frxs vs NVPC:  #####
hhcols <- c("seqnames", "profile_accession", "#match/tLen", "ali_len", "mismatch", "GapOpen", "start", "end", "p1", "p2", "evalue", "score")
HHRepsVsNVPC <- fread("./NVPC/SingleTypeProfs/HHsearch_NVPC_Reps6frx.tsv",sep = "\t",col.names = hhcols,nThread = 10)
HHRepsVsNVPC = HHRepsVsNVPC[wh(HHRepsVsNVPC$evalue < 0.01),]
HHRepsVsNVPC = HHRepsVsNVPC[score >= 8,]
HHRepsVsNVPC$PlenCov <- HHRepsVsNVPC$p2 - HHRepsVsNVPC$p1
w1 <- intersect(wh(HHRepsVsNVPC$PlenCov < 600),wh(HHRepsVsNVPC$PlenCov > 60))
HHRepsVsNVPC <- HHRepsVsNVPC[w1,]
HHRepsVsNVPC <- cbind(HHRepsVsNVPC,str_split_fixed(string = HHRepsVsNVPC$seqnames,pattern = fixed("."),n = 2))
HHRepsVsNVPC <- Rename1Col(HHRepsVsNVPC,"V1","ND")
HHRepsVsNVPC <- Rename1Col(HHRepsVsNVPC,"V2","frame")
HHRepsVsNVPC$frame <- as.integer(HHRepsVsNVPC$frame)
HHRepsVsNVPC$strand <- "+"
HHRepsVsNVPC$strand[HHRepsVsNVPC$frame < 0 ] <- "-"

HHRepsVsNVPC$Analysis_Type <- "NVPC"
HHRepsVsNVPC$Classified <- 3
HHRepsVsNVPC$PolyPort_Problematic = F
HHRepsVsNVPC$Comment <- "Profiles made from reseeded matches."
HHRepsVsNVPC <- merge(HHRepsVsNVPC,NVPC_Profs,by="profile_accession",all.x=T,all.y=F)
HHRepsVsNVPC <- Rename1Col(HHRepsVsNVPC,"q1","start")
HHRepsVsNVPC <- Rename1Col(HHRepsVsNVPC,"q2","end")

framedf <- HHRepsVsNVPC

setorderv(framedf,cols  =  c("seqnames","Classified","score","evalue"),order = c(-1,-1,-1,1))
framedf_asgr <- plyranges::as_granges(framedf) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 50, framedf_asgr,ignore.strand = T,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}


HHRepsVsNVPCDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(framedf_asgr,stringsAsFactors = F)),factors.as.char  =  T))
# Manual inspection - something is off with these results. The profile-query pairs here tend to come from longer matches/alignments, compared to the same pair of query-profile in the hmmsearch. Spot checking a few, seems the hhsearch longer matches are over extended.
# TL;DR - the HHsearch matches won't be intigrated into the bulk matches for now. 



########## Consolidate Pfam34A, RnaVirDB200205, ECOD, CATH, CDD ...  ########## 
AllTMP = rbind(PfamDF, EcodDF, fill = T)
AllTMP = rbind(AllTMP, Cathdf,fill = T)
AllTMP = rbind(AllTMP, InterDomsdf, fill = T)
AllTMP = rbind(AllTMP, PRINTsdf, fill = T)
AllTMP = rbind(AllTMP, CDDF, fill = T)
AllTMP = rbind(AllTMP, SCOPeDF, fill = T)
AllTMP = rbind(AllTMP, SPDF, fill = T)
AllTMP = rbind(AllTMP, LysDF, fill = T)
AllTMP = rbind(AllTMP, p0002MA, fill = T)
AllTMP = distinct(rbind(AllTMP, Rna200205df, fill = T))
AllTMP = distinct(`dropcols<-`(AllTMP,c("PolyPort_Problematic","DB","width","Classified","Analysis_Type","cath.superfamily","profile","IPR_accession", "IPR_annotations","name", "profile_name")))
AllTMP = distinct(`dropcols<-`(AllTMP,c("strand",setdiff(SharedCols(CM12,AllTMP),"profile_accession"),"width")))

AllTMP = merge(AllTMP,data.frame(CM12),by="profile_accession",all.x=T,all.Y=F)
which.na(AllTMP$Classified)
AllTMP$Classified = as.num(AllTMP$Classified)# AllTMP$Dwhich.na(AllTMP$Classified)
AllTMP$Classified[grep(pattern = "lys",x = AllTMP$matchID,ignore.case = T)] = 2
gc()
setDT(AllTMP)

setorderv(AllTMP,cols  =  c("seqnames","Classified","score"),order = c(-1,-1,-1))

framedf_asgr <- plyranges::as_granges(AllTMP) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()
for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = gpp, minoverlap = 50, framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1) != 0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 90,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 30, framedf_asgr,ignore.strand = F,type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

gc()
AllmsDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(framedf_asgr,stringsAsFactors = F)),factors.as.char  =  T))
AllmsDF = distinct(`dropcols<-`(AllmsDF,c(setdiff(SharedCols(CM12,AllmsDF),c("profile_accession")),"width")))
AllmsDF = Rename1Col(AllmsDF,"seqnames","ORFID")

AllDF = distinct(`dropcols<-`(AllmsDF,c(setdiff(SharedCols(CM12,AllmsDF),c("profile_accession")),"width")))

AllDF <- Rename1Col(AllDF,"start","ORF_q1")
AllDF <- Rename1Col(AllDF,"end","ORF_q2")
AllDF <- Rename1Col(AllDF,"qL","ORF_qL")
AllDF$strand <- NULL

AllDF2 <- merge(AllDF,ORFsInfo[,c("start","end","ORFID","strand","seqid")],by="ORFID",all.x=T,all.y=F)
AllDF2 <- Rename1Col(AllDF2,"seqid","ND")
AllDF2 <- merge(AllDF2,IDFT[,c("Length","ND")],by="ND",all.x=T,all.y=F)
AllDF2 <- Rename1Col(AllDF2,"Length","Contig_Length")

AllDF2 <- Rename1Col(AllDF2,"start","Nuc_Start")
AllDF2 <- Rename1Col(AllDF2,"end","Nuc_End")
fwrite(sep = '\t',file = "AllDF.23112021.tsv",x = AllDF2,quote = T,verbose = T)

fwrite(sep = '\t',file = "LysDF2.23112021.tsv",x = LysDF2,quote = T,verbose = T)


##
AllDF = merge(AllDF,CM12,by="profile_accession",all.x=T,all.y=F)
AllDF$freq = NULL
AllDF$Analysis_Type[AllDF$profile_accession %in% LysDF$profile_accession] = "Lysis-LysDB"
w1 = intersect(x = which.na(AllDF$New_Name),y = grep(pattern = "Lys",x = AllDF$matchID,fixed = T))
AllDF$New_Name[w1] = "Lysis?"

View(AllDF[which.na(AllDF$Analysis_Type),])

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "CDD",x = AllDF$matchID,fixed = T))
AllDF$Analysis_Type[w1] = "CDD3.19"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "CATH",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "CATH"
AllDF$Analysis_Type[which(AllDF$Analysis_Type == "Cath")] = "CATH"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "IntDom",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "InterDomains"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "RVD25",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "RNAVirDB200205"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "ECOD",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "ECOD"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "SCOPe",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "SCOPe"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "P34",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "Pfam34A"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "Splitted_Prof",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "Splitted_Profs"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "PRINTS",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "PRINTS"

TmpAllDF <- AllDF

AllDF <- filter(AllDF,!(ORFID %in% NVPCDF$ORFID))
AllDF <- rbind.fill(AllDF,NVPCDF)

AllDF <- Rename1Col(AllDF, "ORFID", "seqnames")
AllDF$Classified[wh(AllDF$Classified > 0 )] = 1 
AllDF$Classified[wh(AllDF$Classified <= 0 )] = -1 
setorderv(AllDF,cols  =  c("seqnames","Classified","score"),order = c(-1,-1,-1))

# AllDF <- AllDF[-wh(AllDF$evalue>0.01),]
framedf_asgr <- plyranges::as_granges(AllDF) #use namespace, avoid silliness.
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F, type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()

TMP1 = (GenomicRanges::findOverlaps(framedf_asgr,ignore.strand = F,type ='equal', select  =  "first",drop.self = F))
w1 = which.na(TMP1)
if(length(w1)!=0){TMP1[w1] = w1}
framedf_asgr = framedf_asgr[unique(TMP1)]
gc()

for(gpp in c(-1,10,93)){
  TMP1 = (GenomicRanges::findOverlaps(maxgap = gpp, minoverlap = 120, framedf_asgr, ignore.strand = F, type ='equal', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:3)){
  TMP1 = (GenomicRanges::findOverlaps(framedf_asgr, ignore.strand = F, type ='within', select  =  "first", drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}

for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps( maxgap = 150,minoverlap = 30, framedf_asgr,ignore.strand = F,type ='within', select  =  "first",drop.self = T))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}


for(i in c(1:2)){
  TMP1 = (GenomicRanges::findOverlaps(minoverlap = 250, framedf_asgr, ignore.strand = F, type ='any', select  =  "first",drop.self = F))
  w1 = which.na(TMP1)
  if(length(w1)!=0){TMP1[w1] = w1}
  framedf_asgr = framedf_asgr[unique(TMP1)]
  gc()
}




gc()
AllmsDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(framedf_asgr,stringsAsFactors = F)),factors.as.char  =  T))
AllmsDF = distinct(`dropcols<-`(AllmsDF,c(setdiff(SharedCols(CM12,AllmsDF),c("profile_accession")),"width")))

AllmsDF = Rename1Col(AllmsDF,"seqnames","ORFID")
AllDF = distinct(`dropcols<-`(AllmsDF,c(setdiff(SharedCols(CM12,AllmsDF),c("profile_accession")),"width")))

AllDF = merge(AllDF,CM12,by="profile_accession",all.x=T,all.y=F)
AllDF$freq = NULL
AllDF$Analysis_Type[AllDF$profile_accession %in% na.omit(LysDF$profile_accession)] = "Lysis-LysDB"
w1 = intersect(x = which.na(AllDF$New_Name),y = grep(pattern = "Lys",x = AllDF$matchID,fixed = T))
AllDF$New_Name[w1] = "Lysis?"

# View(AllDF[which.na(AllDF$Analysis_Type),])

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "CDD",x = AllDF$matchID,fixed = T))
AllDF$Analysis_Type[w1] = "CDD3.19"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "CATH",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "CATH"
AllDF$Analysis_Type[which(AllDF$Analysis_Type == "Cath")] = "CATH"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "IntDom",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "InterDomains"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "RVD25",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "RNAVirDB200205"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "ECOD",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "ECOD"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "SCOPe",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "SCOPe"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "P34",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "Pfam34A"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "Splitted_Prof",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "Splitted_Profs"

w1 = intersect(x = which.na(AllDF$Analysis_Type),y = grep(pattern = "PRINTS",x = AllDF$matchID,fixed = F))
AllDF$Analysis_Type[w1] = "PRINTS"


AllDF$DB <- NULL
AllDF$Ver <- NULL

oldDF <- fread(na.strings = "",input = "AllDF.15112021.tsv")
oldDF <- filter(oldDF,ORFID %in% c(NScodexNVPC$ORFID,setdiff(oldDF$ORFID,AllDF$ORFID)))
AllDF <- filter(AllDF,!(ORFID %in% NScodexNVPC$ORFID))
setdiff(oldDF$ORFID,AllDF$ORFID)

tmpcols <- SharedCols(oldDF,AllDF)
AllDF <- rbind.fill(data.frame(AllDF)[,tmpcols],data.frame(oldDF)[,tmpcols])
saveRDS(object = AllDF,"AllDF.22112021.RDS")
# AllDF <- readRDS("AllDF.22112021.RDS")
fwrite(sep = '\t',file = "AllDF.22112021.tsv",x = AllDF,quote = T,verbose = T)
colnames(RdRpDF) <- c("ND","ORFID","start","end", "profile_accession", "evalue", "score", "p1", "p2", "pL", "strand")
RdRpDF$New_Name <- "RdRp"
RdRpDF$Analysis_Type <- "RdRpDB"
RdRpDF$Classified <- 1 

OrfsInfo$strand
AllDF <- rbind.fill(data.frame(AllDF),data.frame(RdRpDF))
AllDF$ali_Qcov_len <- AllDF$end - AllDF$start

w1 <- wh(AllDF$ORFID %in% RdRpDF$ORFID)
AllDF <- AllDF[-w1,]
w1 <- which.max(AllDF$ali_Qcov_len)
View(AllDF[w1,])
pdf("AllDF_ali_Qcov_len_Hist.pdf", 
    width = 10, 
    height = 10)  
ggplot(data.frame("size" = AllDF$ali_Qcov_len),aes(size)) +  geom_histogram(bins = 1500) + scale_y_log10()
dev.off()
w1 <- wh(AllDF$ali_Qcov_len < 20)
AllDF <- AllDF[-w1,]
AllDF[which.min(na.omit(AllDF$ali_Qcov_len)),]

{ # 22.11.2021 - Additional annotations from Mart.
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("AGQ47773-1-741","AQU42764-9-750","APG79334-949-1921","YP_009052470-201-645"))] <- "RdRp"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("ANW72252-660-846","YP_004598981-675-879","ANW72253-663-850"))] <- "Torsin-1A-interacting protein 1"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("ANW72253-1339-1651","ANW72252-1339-1651","YP_004598981-1387-1700"))] <- "PRO−Chy"
  AllDF$New_Name[wh(AllDF$New_Name %in% c("Cap_MTase_GTase"))] <- "Cap_MTase−GTase"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("AJZ68872-1-212"))] <- "Cap_EndoN"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("AEX65765-21-494"))] <- "Env_CIF"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_002308458-40-578","YP_002905338-6-577","NP_049359-603-1101"))] <- "Env_CIIIF"
   
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_001661452-6000-6255"))] <- "CoV_NSP14_ExoN"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("NP_690824-1-331"))] <- "Helicase_SF4"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("NP_524563-1-792"))] <- "CP_Cysto_T1"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("AEX65762-1-524"))] <- "CP_Borna"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_002332933-6-237"))] <- "CP_Flexi/Phlebo"
  
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("NP_041090-1-243"))] <- "PRO-Pap/vOTU"
  AllDF$Comment[wh(AllDF$profile_accession %in% c("NP_041090-1-243"))] <- "Cys protease (PF01830.21; MEROPS family C7) distantly related to papain-like proteases"
  
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("NP_620651-109-308"))] <- "Allexi_40kDa"

  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_009342047-252-391"))] <- "Macro"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_009252304-1900-2087"))] <- "Vpg_Poty"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_009344970-1996-2231","56112__d.144.1__d1tqia2","YP_009344970-1621-1859"))] <- "STKc"
  AllDF$Comment[wh(AllDF$profile_accession %in% c("YP_009344970-1996-2231","56112__d.144.1__d1tqia2","YP_009344970-1621-1859"))] <- "Serine/Threonine Kinase, cAMP-dependent protein kinase"

  AllDF$New_Name[wh(AllDF$profile_accession %in% c("YP_009094172-31-568"))] <- "Env_HN"
  AllDF$Comment[wh(AllDF$profile_accession %in% c("YP_009094172-31-568"))] <- "Envelope protein, hemagglutinin-neuraminidase"
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("PF01806.19"))] <- "Paramyxo_P"
  AllDF$Comment[wh(AllDF$profile_accession %in% c("PF01806.19"))] <- "phosphoprotein P of paramyxoviruses"
  
  
  AllDF$New_Name[wh(AllDF$profile_accession %in% c("COG0741","cd00736","COG4678","3vwoA01","4mphA00","2vo9A01","1lbuA02","PF05838"))] <- "LYZF1"

  AllDF$New_Name[wh(AllDF$New_Name %in% c("Lysis_Levi"))] <- "Lysis_sgl"
  AllDF$New_Name[wh(AllDF$New_Name %in% c("PRO-M34 (lysis?)"))] <- "Lysin_M34"
  AllDF$New_Name[wh(AllDF$New_Name %in% c("PRO-M15 (lysis?)"))] <- "Lysin_M15"
  
  AllDF$New_Name[wh(AllDF$New_Name %in% c("Lysis_M23 family peptidase","Lysis_PRO-M23"))] <- "Lysin_M23"
  
  
  AllDF$Classified[wh(AllDF$profile_accession %in% c("YP_009344970-1621-1859","NP_524563-1-792","NP_690824-1-331","YP_009344970-1996-2231","YP_001661452-6000-6255","NP_049359-603-1101","YP_002308458-40-578","YP_002905338-6-577","YP_009252304-1900-2087","YP_009052470-201-645","YP_009342047-252-391","APG79334-949-1921","AEX65765-21-494","AJZ68872-1-212","ANW72253-1339-1651","ANW72252-1339-1651","YP_004598981-1387-1700","AGQ47773-1-741","ANW72252-660-846","AQU42764-9-750","YP_004598981-675-879","ANW72253-663-850"))] <- 3
  AllDF$PolyPort_Problematic[AllDF$profile_accession %in% PolyProts] = T
  AllDF <- filter(AllDF,profile_accession !="cluster_003204")
  
}
CM13 <- distinct(AllDF[,SharedCols(CM12,AllDF)])
WriteXlsx(CM13,"CM13.xlsx")
# TmpAllDF <- AllDF
# AllDF <- TmpAllDF
# AllDF$ORFID %in% 


# tmmp = distinct(data.frame(AllDF)[,SharedCols(AllDF,CM11)])
# tmmp = distinct(data.frame(AllDF)[,SharedCols(AllDF,CM11)])
# tmmp$profile_name = tmmp$profile_accession
# w1 = which(tmmp$Analysis_Type == "SCOPe")
# tmmp$profile_name[w1] = paste(tmmp$SCOPe_Function_ID[w1],tmmp$SCOPe_Class[w1],tmmp$SCOPe_PDBID[w1],sep = "__")
# tmmp$profile_accession[w1] = tmmp$profile_name[w1]
saveRDS(object = AllDF,"AllDF.23112021.RDS")

# WriteXlsx(dfdt = AllDF,"AllDF.05102021.xlsx")
fwrite(sep = '\t',file = "AllDF.23112021.tsv",x = AllDF,quote = T,verbose = T)

rm(AllTMP,tmmp,w1,framedf_asgr,framedf)
gc()
remcols = c( "eff_nseq", "relent","info","p_relE","compKL","Analysis_Type.y","freq","Analysis_Type")
TMPHMMstats15 = `dropcols<-`(HMMstats15,remcols)
PfamDF2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(PfamDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
EcodDF2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(EcodDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
Cathdf2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(Cathdf,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
InterDomsdf2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(InterDomsdf,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
PRINTsdf2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(PRINTsdf,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
CDDF2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(CDDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
SCOPeDF2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(SCOPeDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))
Rna200205df2 = Rename1Col(newcolnm = "ORFID",colnm = "seqnames",dfdt = RemovesSparseols(merge(Rna200205df,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)))


PhobTMModdf = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(TMP_gr,stringsAsFactors = F)),factors.as.char  =  T))

PhobTMModdf$evalue =NULL
PhobTMModdf$matchID =NULL
PhobTMModdf$DB =NULL
PhobTMModdf$IPR_accession =NULL
PhobTMModdf$evalue =NULL
PhobTMModdf$ND = PhobTMModdf$seqnames
PhobTMModdf$seqnames = NULL
PhobDF = filter(PhobTMModdf,Analysis_Type =="Phobius")
MobiDBDF = filter(PhobTMModdf,Analysis_Type =="MobiDBLite")
TMHMMDF = filter(PhobTMModdf,Analysis_Type =="TMHMM")

EcodDF2 = merge(EcodDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
Cathdf2 = merge(Cathdf,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
InterDomsdf2 = merge(InterDomsdf,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
PRINTsdf2 = merge(PRINTsdf,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
CDDF2 = merge(CDDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
SCOPeDF2 = merge(SCOPeDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
SPDF2 = merge(SPDF,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)
Rna200205df2 = merge(Rna200205df,TMPHMMstats15,by="profile_accession",all.x=T,all.y=F)

WriteXlsx(EcodDF2,"ECOD_Hitbl.xlsx")
WriteXlsx(Cathdf2,"CATH_Hitbl.xlsx")
WriteXlsx(InterDomsdf2,"InterDoms_Hitbl.xlsx")
WriteXlsx(PRINTsdf2,"PRINTs_Hitbl.xlsx")
WriteXlsx(CDDF2,"CDD_Hitbl.xlsx")
WriteXlsx(SCOPeDF2,"SCOPe_Hitbl.xlsx")
WriteXlsx(SPDF2,"SPD_Hitbl.xlsx")
WriteXlsx(Rna200205df2,"RNA200205_Hitbl.xlsx")
write_excel_csv2(num_threads = 4,x = Rna200205df2,file = "RNA200205_Hitbl.tsv",delim = '\t',na = "",col_names = T,progress = T)
WriteXlsx(PhobDF,"Phobius_Hitbl.xlsx")
WriteXlsx(MobiDBDF,"MobiDBLite_Hitbl.xlsx")
WriteXlsx(TMHMMDF,"TMHMM_Hitbl.xlsx")

WriteXlsx(AllDF,"AllDF.05102021.tsv")

fwrite(sep = '\t',file = "AllDF.05102021.2.tsv",x = AllDF,quote = T,verbose = T)
WriteXlsx(dfdt = AllTbl,filepath = "Alltbl.0310.2021.xlsx")
WriteXlsx(dfdt = CM12,filepath = "CM11.xlsx")


ORFsInfo_gr$Analysis_Type = "ORF_Prediction"
AllGR = c(ORFsInfo_gr,Allms_gr)
testttt = unique(union(as.char(ORFsInfo_gr@seqnames),as.char(Allms_gr@seqnames)))
setdiff(IDFT$Contig,testttt)
setdiff(testttt, IDFT$ND)

AllDFxIDFT = merge(AllDF,Rename1Col(IDFT,"Contig","ND"),by="ND",all.x=T,all.y=F)
AllDFxIDFT = filter(AllDFxIDFT, evalue<0.001)
TMP2 = plyr::count(AllDFxIDFT,vars = c("Phylum","New_Name"))
TMP2 = plyr::count(distinct(data.frame(AllDFxIDFT)[,c("New_Name","ND","profile_accession")]),vars = c("profile_accession"))
CM12$freq.y = NULL
CM12 = merge(TMP2[,c("freq","profile_accession")],CM12,by="profile_accession",all.x=T,all.y=F)
CM12$freq[which.na(CM12$freq)] = 0 
WriteXlsx(dfdt = CM12,filepath = "CM13.xlsx")

 # 1. Extract all uncovered regions

# 2. Add TMHMM, PHOBIUS and Moddblite.
colnames(PhobTMModdf) = c("ORFID", "Analysis_Type", "profile_accession", "profile", "start", "end", "IPR_accession","IPR_annotations","evalue","ali_Qcov_len", "seqnames","strand","matchID","DB" )
TMP_gr = plyranges::as_granges(PhobTMModdf)
AllGR = c(AllGR,TMP_gr)
AllDF = setDT(BBmisc::convertDataFrameCols((GenomicTuples::as.data.frame(AllGR,stringsAsFactors = F)),factors.as.char  =  T))
AllDF = Rename1Col(AllDF,"seqnames","ND")
WriteXlsx(dfdt = AllDF,filepath = "AllDF.03102021.xlsx")
fwrite(sep = '\t',file = "AllDF.03102021.tsv",x = AllTbl,quote = T,verbose = T)
saveRDS(AllDF,"AllDF.03102021.RDS")
# saveRDS(CM12,"CM12.RDS")

WriteWolfTbl(AllDF,"AllDF.03102021.tsv")
WriteWolfTbl(TMP2,"DB_hit_count.tsv")
WriteWolfTbl(TMP1,"FunClass_hit_count.tsv")
WriteWolfTbl(TMP3,"PhylumXFunClass_hit_count.tsv")
WriteWolfTbl(TMP4,"PhylumXProfileName_hit_count.tsv")
WriteWolfTbl(TMP5,"PhylumXProfileID_hit_count.tsv")

WriteWolfTbl(filter(AllmsDFxIDFT,ND %in% filter(AllDF,FunClass == "Lysis")$ND) ,"LysDF.tsv")
writeXStringSet(IDFT.fasta[unique(filter(AllDF,FunClass == "Lysis")$ND)],"LysFa.fasta",width=20001)

AllGRlist = plyranges::group_by(.data = AllGR, ORFID)
tmp2 = plyranges::compute_coverage(AllGR)
tmp2 = GenomicRanges::grglist(AllGR)

AllGRlist = GRangesList(AllGR)
AllGFF = rtracklayer::asGFF(AllGRlist)
write_gff3(x = AllGR,file = "AllGR.gff2")

AllmsDFxIDFT = merge(IDFT,Rename1Col(AllmsDF,"seqnames","ND"),by="ND",all.x=T,all.y=F)
AllmsDFxIDFT = `dropcols<-`(AllmsDFxIDFT,c("Note","pCoverage", "evalue" ,"DB","Analysis_Type" ,"end","strand","p1","p2" , "matchID", "Classified","nuc.cls9090","cls9595_id","RdRp_ID","UVA","nseq","M", "eff_nseq", "relent", "info", "p_relE", "compKL"))
countDF = plyr::count(AllmsDFxIDFT[,c("Order","profile","New_Name")])
# WriteWolfTbl(AllmsDFxIDFT,"AllmsDFxIDFT.tsv")
# WriteWolfTbl(Rename1Col(AllmsDF,"seqnames","ND"),"AllmsDF.tsv")

countDF 
TMP4 = plyr::count(AllmsDFxIDFT[,c("Phylum","profile")])
TMP5 = plyr::count(AllmsDFxIDFT[,c("Phylum","profile_accession")])

rm(TMP1,TMP2,framedf,framedf_asgr,PfDBms,RDBms,HairEyeColor,HPRDBms,HPVA,HRVA,w2,w4,w3,w1)

####### Profile Dereplication ####### 
Profs2709 = distinct(readxl::read_xlsx("/home/neri/Downloads/ProfCounts0209tax_20210917.xlsx",na = "NA",col_types = "text"))
Profs2709$...6 = NULL
Profs2709$...5 = NULL

w1 =-which.na(HMMstats15$SCOPe_Class)
HMMstats16 = HMMstats15

HMMstats16$profile_accession[w1] = paste(HMMstats16$SCOPe_Function_ID[w1],HMMstats16$SCOPe_Class[w1],HMMstats16$SCOPe_PDBID[w1],sep = "__")

HMMstats16$New_Name[which(HMMstats16$New_Name == "-")] = NA

tmpdfff5 = merge(x = Profs2709, y =  HMMstats16, by = "profile_accession",all.x=T, all.y=F)
compdf = tmpdfff5
compdf[,c("Comp_New_Name","Comp_Comment")] = F
setDF(compdf)
for (colnym in c("New_Name","Comment")) {
  compdf[,p0("Comp_",colnym)] = (compdf[,p0(colnym,".x")] == compdf[,p0(colnym,".y")])
}

compdf$Comp_Comment[intersect(which.na(compdf$Comment.x),which.na(compdf$Comment.y))] = T
compdf$Comp_New_Name[intersect(which.na(compdf$New_Name.x),which.na(compdf$New_Name.y))] = T

View(compdf[which.na(compdf$Comp_Comment),c("profile_accession","New_Name.x","Comment.x","New_Name.y","Comment.y","Comp_New_Name","Comp_Comment")])

setdiff(unique(CM12$New_Name),compdf$New_Name)
WriteXlsx(compdf,"compdf.xlsx")#
compdf1 = ReadXlsx(filepath = "compdf.xlsx")

# compdf2 = compdf
# compdf2$Tsss = apply(X = compdf,MARGIN = 1,FUN = function(x) (na.omit(unique(x["New_Name.x"],x["New_Name.y"]))))
# compdf2$Tsss = as.char(compdf2$Tsss)
# compdf2$Tsss[which(compdf2$Tsss == "character(0)")] = NA
# View(compdf2[,c("profile_accession","Tsss","New_Name.x","Comment.x","New_Name.y","Comment.y","Comp_New_Name","Comp_Comment")])
# compdf2$New_Name = compdf2$Tsss
# compdf2$Tsss = apply(X = compdf,MARGIN = 1,FUN = function(x) (na.omit(unique(x["Comment.x"],x["Comment.y"]))))
# compdf2$Tsss = as.char(compdf2$Tsss)
# compdf2$Tsss[which(compdf2$Tsss == "character(0)")] = NA
# View(compdf2[,c("profile_accession","Tsss","New_Name.x","Comment.x","New_Name.y","Comment.y","Comp_New_Name","Comp_Comment")])
# compdf2$Comment = compdf2$Tsss
# 
# # 
# 
# 
# compdf3 = `dropcols<-`(compdf2,c("Tsss","New_Name.x","Comment.x","New_Name.y","Comment.y","Comp_New_Name","Comp_Comment"))
compdf3 = compdf1
which.na(compdf3$profile_name)
compdf4 = merge(compdf3,Profs2208_mod[,c("profile_accession","New_Name")],by="profile_accession",all.x=T,all.y=F)
compdf4$profile_name[which.na(compdf4$profile_name)] = compdf4$name.x[which.na(compdf4$profile_name)]

compdf4[which(compdf4$name.y != compdf4$name.y),]
compdf4 = `dropcols<-`(compdf4,c("name.y","name.x"))


sharedcols = intersect(colnames(compdf4),colnames(Profs2208_mod))
sharedcols = setdiff(sharedcols,c("profile_accession","Comment","New_Name","name"))

compdf5 = merge(x = compdf4, y =  Profs2208_mod[,c("profile_accession",sharedcols)], by = "profile_accession",all.x=T, all.y=F)
compdf5[,p0("Comp_",sharedcols)] = F

setDF(compdf5)
compdf5[,(sharedcols)] = NA

for (colnym in sharedcols) {
  compdf5[,p0("Comp_",colnym)] = (compdf5[,p0(colnym,".x")] == compdf5[,p0(colnym,".y")])
  nass = intersect(which.na(compdf5[,p0(colnym,".x")]),which.na(compdf5[,p0(colnym,".y")]))
  compdf5[nass, p0("Comp_",colnym)] = T 
  nass = which.na(compdf5[,p0("Comp_",colnym)])
  if(length(nass)!=0){
  compdf5[nass, p0("Comp_",colnym)] = F 
  }
  
  compdf5[,colnym] = compdf5[,p0(colnym,".x")]
  w1 = which.na(compdf5[,p0(colnym,".x")])
  w2 = which.Not.na(compdf5[,p0(colnym,".y")])
  w3 = intersect(w2,w1)
  if(length(w3)!=0){
  compdf5[w3,colnym] = compdf5[w3,p0(colnym,".y")]
  }
}

compdf5[,"AllCompsTrue"] = apply(X = compdf5[,p0("Comp_",sharedcols)],MARGIN = 1,AllTrue)
cols2drop = c(p0("Comp_",sharedcols),c(p0(sharedcols,".y"),(p0(sharedcols,".x"))))
# WonderDF = `dropcols<-`(compdf5[which(compdf5[,"AllCompsTrue"]),],cols2drop)
WonderDF = `dropcols<-`(compdf5,cols2drop)
WonderDF = distinct(WonderDF)
WonderDF[which.duplicated(WonderDF$profile_accession),]
data.table::fwrite(WonderDF,"WonderDF.tsv",quote = T)
WonderDF = readxl::read_xlsx("/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/WonderDF.xlsx",na = "NA",col_types = "text")
# 
# compdf5 = merge(compdf4,WonderDF,by="profile_accession",all.x=T,all.y=F)
WonderDF = compdf4
WonderDF$Ver = "03.10.2021"
# WonderDF$New_Name = WonderDF$New_Name.x
# WonderDF$New_Name.x = NULL
# WonderDF$Analysis_Type = WonderDF$Analysis_Type.x
# which.na()

{
  {
    SCOPe_Desc$profile_accession
  }
  
  HMMstats15$Analysis_Type[which.Not.na(HMMstats15$SCOPe_PDBID)] = "SCOPe"
  w1 = which(HMMstats15$Analysis_Type == "SCOPe")
  HMMstats15$name1 = NA
  HMMstats15$name1[w1] = HMMstats15$profile_accession[w1]
  HMMstats15$profile_accession[w1] = HMMstats15$name[w1]
  WonderDF = merge(WonderDF,HMMstats15[c("profile_accession","name1")],all.x=T,all.y=F)
  FatWonderDF = gtools::smartbind(data.frame(WonderDF),data.frame(HMMstats15[which(!(HMMstats15$profile_accession %in% WonderDF$profile_accession)),]))
}

# FatWonderDF = rbind.fill(data.frame(WonderDF),data.frame(Profs2208[which(!(Profs2208$profile_accession %in% WonderDF$profile_accession)),]))
which.duplicated(FatWonderDF$profile_accession)
FatWonderDF[which.duplicated(FatWonderDF$profile_accession)]

FatWonderDF = rbind.fill(data.frame(FatWonderDF),data.frame(HMMstats15[which(!(HMMstats15$profile_accession %in% WonderDF$profile_accession)),]))
FatWonderDF[which.duplicated(FatWonderDF$profile_accession),]
FatWonderDF = `droprows<-`(FatWonderDF,82816)
# clan_membership = scan(what = "\n",sep = "\n",file = "/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/clan_membership.txt")
# clan_membership = data.frame("mems" = clan_membership,
#                              "ClanID"  = p0("Clan_",pad(str = c(1:length(clan_membership)),pad ="0",side = "left",use_length = T,width = nchar(length(clan_membership)))) )
# 
# clan_membership = TabUnroll(dfdt = clan_membership, sep = ",", colnym = "mems",NewColNym = "profile_accession")
# clan_membership = `dropcols<-`(clan_membership,"mems")
# 
# setdiff(clan_membership$profile_accession,HMMstats15$profile_accession)
# clan_membership = merge(clan_membership,FatWonderDF[,c("profile_accession","New_Name")],by="profile_accession",all.x=T,all.y=F)
# # 
# # clan_membership_cout = plyr::count(clan_membership,vars = c("ClanID","New_Name"))
# # 
# # UniqClans = clan_membership2[-which.duplicated(clan_membership2$ClanID),]
# # DupedClans = clan_membership2[which.duplicated(clan_membership2$ClanID),]
# # # setorderv(setDT(DDllyst),c("ND.frame","score","pident","evalue"),order = c(-1,-1,-1,1))
# # # DDllyst = unique(DDllyst, incomparables=FALSE, fromLast=FALSE,by="ND.frame")
# # 
# 
# w1 = intersect(which.duplicated(clan_membership2$ClanID,returnalldups = T),which.na(clan_membership2$New_Name))
# 
# View(clan_membership2[w1,])
# 
# clan_membership2 = clan_membership %>% group_by(ClanID) %>% summarise(
#   MemsList= list(unlist(list(unique(profile_accession)))),
#   MemsChar= toString(unlist(MemsList)),
#   Nmems = lengths(MemsList),
#   NameList= list(unlist(list(unique(na.omit(New_Name))))),
#   NamesChar= toString(unlist(NameList)),
#   Nnames = lengths(NameList)#,
#   # New_Name = New_Name
# )
# clan_membership3 = `dropcols<-`(clan_membership2,c("MemsList","NameList"))
# clan_membership3 = TabUnroll(dfdt = clan_membership3,sep = ", ",colnym = "MemsChar",NewColNym = "profile_accession")
# clan_membership3 = merge(clan_membership3,clan_membership[,c("New_Name","profile_accession")],by = "profile_accession", all.x=T,all.y=F)
# 
# xlsx::write.xlsx(clan_membership3,"clan_membership3.xlsx")# 
# clan_membership_3mod = readxl::read_xlsx("/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/clan_membership3_mod.xlsx",na = "NA",col_types = "text")
# clan_membership_3mod = `dropcols<-`(clan_membership_3mod,c("Nmems","Nnames","NamesChar"))
# TempWonderDF = merge(WonderDF,clan_membership_3mod,by = "profile_accession",all.x=T,all.y=F)
# 
# sharedcols = c("New_Name")
# 
# 
# TempWonderDF[,p0("Comp_",sharedcols)] = F
# 
# setDF(TempWonderDF)
# TempWonderDF[,(sharedcols)] = NA
# 
# for (colnym in sharedcols) {
#   TempWonderDF[,p0("Comp_",colnym)] = (TempWonderDF[,p0(colnym,".x")] == TempWonderDF[,p0(colnym,".y")])
#   nass = intersect(which.na(TempWonderDF[,p0(colnym,".x")]),which.na(TempWonderDF[,p0(colnym,".y")]))
#   TempWonderDF[nass, p0("Comp_",colnym)] = T 
#   nass = which.na(TempWonderDF[,p0("Comp_",colnym)])
#   if(length(nass)!=0){
#     TempWonderDF[nass, p0("Comp_",colnym)] = F 
#   }
#   
#   TempWonderDF[,colnym] = TempWonderDF[,p0(colnym,".x")]
#   w1 = which.na(TempWonderDF[,p0(colnym,".x")])
#   w2 = which.Not.na(TempWonderDF[,p0(colnym,".y")])
#   w3 = intersect(w2,w1)
#   if(length(w3)!=0){
#     TempWonderDF[w3,colnym] = TempWonderDF[w3,p0(colnym,".y")]
#   }
# }
# 
# View(TempWonderDF[,c(sharedcols,"New_Name.x","New_Name.y")])
# 
# 
# 
# clan_membership5 = TempWonderDF %>% group_by(ClanID) %>% summarise(
#   MemsList= list(unlist(list(unique(profile_accession)))),
#   MemsChar= toString(unlist(MemsList)),
#   Nmems = lengths(MemsList),
#   NameList= list(unlist(list(unique(na.omit(New_Name))))),
#   NamesChar= toString(unlist(NameList)),
#   Nnames = lengths(NameList)#,
#   # New_Name = New_Name
# )
# clan_membership6 = `dropcols<-`(clan_membership5,c("MemsList","NameList"))
# clan_membership6 = TabUnroll(dfdt = clan_membership6,sep = ", ",colnym = "MemsChar",NewColNym = "profile_accession")
# clan_membership6 = merge(clan_membership6,TempWonderDF[,c("New_Name","profile_accession")],by = "profile_accession", all.x=T,all.y=F)
# clan_membership6 = clan_membership6[-which.na(clan_membership6$ClanID),]
# xlsx::write.xlsx(clan_membership6,"clan_membership6.xlsx")# 
# xlsx::write.xlsx(TempWonderDF,"TempWonderDF.xlsx")# 
# 
# 
# w1 = which.na(Profs27092$New_Name)
# w2 = intersect(Profs27092$profile_accession[w1], TempWonderDF$profile_accession[which.Not.na(TempWonderDF$New_Name)])
# sum(as.integer(na.omit(HMMstats15$freq[which(HMMstats15$profile_accession %in% w2)])))
#   TempWonderDF
# # 
# # 
# # for (colnym in sharedcols) {
# #   compdf5[,colnym] = compdf5[,p0(colnym,".x")]
# #   w1 = which.na(compdf5[,colnym])
# #   
# #   compdf5[,colnym] = (compdf5[,p0(colnym,".x")] == compdf5[,p0(colnym,".y")])
# #   nass = intersect(which.na(compdf5[,p0(colnym,".x")]),which.na(compdf5[,p0(colnym,".y")]))
# #   compdf5[nass, p0("Comp_",colnym)] = T 
# # }
# # 
# # compdf$Comp_Comment[intersect(which.na(compdf$Comment.x),which.na(compdf$Comment.y))] = T
# # compdf$Comp_New_Name[intersect(which.na(compdf$New_Name.x),which.na(compdf$New_Name.y))] = T
# write_lines(TempWonderDF$profile_accession,"tmp_prof_list.txt")
# fwrite(file = "tmp_prof_list.tbl",TempWonderDF[,c("profile_accession","profile_name","Analysis_Type")],quote = T)
# fwrite(sep = '\t',file = "single_domain_profiles.tsv",TempWonderDF[-which.na(TempWonderDF$New_Name),c("profile_accession","profile_name","Analysis_Type","New_Name","Comment")],quote = T)
# TempWonderDF$M =as.num(TempWonderDF$M)
# 
# # write_lines(x = unique(TempWonderDF$profile_accession),file = "all_profile_accessions.lst")


##### Anew #####
ALLprof_clan_membership = scan(what = "\n",sep = "\n",file = "/media/HDD1/uri/RNA_Vir_MTs/RVMT/Annotation/Python/clan_membership.txt")
ALLprof_clan_membership = data.frame("mems" = ALLprof_clan_membership,
                             "ClanID"  = p0("Clan_",pad(str = c(1:length(ALLprof_clan_membership)),pad ="0",side = "left",use_length = T,width = nchar(length(ALLprof_clan_membership)))) )


ALLprof_clan_membership = TabUnroll(dfdt = ALLprof_clan_membership, sep = ",", colnym = "mems",NewColNym = "profile_accession")
ALLprof_clan_membership = `dropcols<-`(ALLprof_clan_membership,"mems")

tmpdf = FatWonderDF[,c("profile_accession","name1","Analysis_Type")]
tmpdf = tmpdf[which(tmpdf$Analysis_Type == "SCOPe"),]
colnames(tmpdf) = c("tmp_profile_accession","profile_accession","Analysis_Type")

tmpdf2 = merge(ALLprof_clan_membership,tmpdf,by="profile_accession",all.x=T,all.y=F)

tmpdf2$tmp_profile_accession[which.na(tmpdf2$tmp_profile_accession)] = tmpdf2$profile_accession[which.na(tmpdf2$tmp_profile_accession)]
ALLprof_clan_membership = tmpdf2[,c("ClanID","tmp_profile_accession")]
colnames(ALLprof_clan_membership) = c("ClanID","profile_accession")

setdiff(ALLprof_clan_membership$profile_accession,FatWonderDF$profile_accession)
which.duplicated(clan_membership$profile_accession)
intersect(clan_membership$profile_accession,FatWonderDF$name)

ALLprof_clan_membership = ALLprof_clan_membership[(!(ALLprof_clan_membership$profile_accession %in% PolyProts)),]

# Stringent_annotation = fread("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Annotation/Python/single_domain_profiles.tsv")
w1 = union(which.Not.na(FatWonderDF$New_Name),which.Not.na(FatWonderDF$Comment))
Stringent_annotation = distinct(FatWonderDF[w1,])
clan_membership = distinct(merge(ALLprof_clan_membership,Stringent_annotation[,c("profile_accession","New_Name")],by="profile_accession",all.x=T,all.y=F))

clan_membership2 = clan_membership %>% group_by(ClanID) %>% summarise(
  MemsList= list(unlist(list(unique(profile_accession)))),
  MemsChar= toString(unlist(MemsList)),
  Nmems = lengths(MemsList),
  NameList= list(unlist(list(unique(na.omit(New_Name))))),
  NamesChar= toString(unlist(NameList)),
  Nnames = lengths(NameList)#,
  # New_Name = New_Name
)

clan_membership_cout = plyr::count(clan_membership,vars = c("ClanID","New_Name"))

UniqClans = clan_membership2[-which.duplicated(clan_membership$ClanID),]
DupedClans = clan_membership2[which.duplicated(clan_membership$ClanID),]
# setorderv(setDT(DDllyst),c("ND.frame","score","pident","evalue"),order = c(-1,-1,-1,1))
# DDllyst = unique(DDllyst, incomparables=FALSE, fromLast=FALSE,by="ND.frame")


w1 = intersect(which.duplicated(clan_membership2$ClanID,returnalldups = T),which.na(clan_membership2$New_Name))

View(clan_membership2[w1,])

clan_membership3 = `dropcols<-`(clan_membership2,c("MemsList","NameList"))
clan_membership3 = TabUnroll(dfdt = clan_membership3,sep = ", ",colnym = "MemsChar",NewColNym = "profile_accession")
clan_membership3 = merge(clan_membership3,clan_membership[,c("New_Name","profile_accession")],by = "profile_accession", all.x=T,all.y=F)

# w1 = which(clan_membership3$Nnames == 1)
w1 = which(clan_membership3$Nnames != 0)

w2 = which(clan_membership_cout$ClanID %in%clan_membership3$ClanID[w1])
AnnotPrecDF = clan_membership_cout[w2,]#data.frame("ClanID" = ClanID,)
AnnotPrecDF$freq[which.na(AnnotPrecDF$New_Name)] = 0
AnnotPrecDF = AnnotPrecDF[-which.na(AnnotPrecDF$New_Name),]
AnnotPrecDF = ddply(data.table(AnnotPrecDF), .(ClanID), transform, AnnotPrec = freq/sum(freq) * 100)
# clan_membership3 = merge(clan_membership3,AnnotPrecDF[,c("ClanID","AnnotPrec")],by=c("ClanID","New_Name"),all.x=T,all.y=F)
clan_membership3$AnnotPrec[which.na(clan_membership3$AnnotPrec)] = 0

WriteXlsx(clan_membership3,"clan_membership3.v3.xlsx")# 
clan_membership_3mod = readxl::read_xlsx("clan_membership3.v3.xlsx",na = "NA",col_types = "text",)
clan_membership_3mod$AnnotPrec = NULL
clan_membership_3mod$...1 =NULL
clan_membership_3mod$NamesChar = NULL
# clan_membership_3mod = distinct(clan_membership_3mod)
# clan_membership_3mod = merge(clan_membership_3mod,clan_membership2[,c("ClanID","NamesChar")],by="ClanID",all.x=T,all.y=F)
# xlsx::write.xlsx(x = clan_membership_3mod,"clan_membership3.v2.xlsx")# 
# 
# clan_membership_3mod = `dropcols<-`(clan_membership_3mod,c("Nmems","Nnames","NamesChar"))
# xlsx::write.xlsx(x = clan_membership_3mod,"clan_membership_3mod.xlsx")# 
# 
# clan_membership_4mod = readxl::read_xlsx("/media/HDD1/uri/RNA_Vir_MTs/V3/Annotations/clan_membership_3mod.xlsx",na = "NA",col_types = "text",)
clan_membership_5mod = distinct(clan_membership_3mod[,c("ClanID","profile_accession","New_Name")])
# clan_membership_6mod = clan_membership_5mod[-which.na(clan_membership_5mod$New_Name),]
clan_membership_5mod$Ver = "03.10.2021"
clan_membership_5mod = merge(clan_membership_5mod,FatWonderDF[,c("profile_accession","profile_name","Comment","name","SCOPe_Function_ID", "SCOPe_Class","SCOPe_PDBID","SCOPe_Desc","profile", "CDD_name", "CDD_desc", "CATH_name","nseq", "eff_nseq", "M", "relent","info","p_relE","compKL","name1")],by="profile_accession",all.x=T,all.y=F)
# clan_membership_5mod[which.duplicated(clan_membership_5mod$profile_accession),]
# clan_membership_5 = clan_membership_5mod[-4480,]
clan_membership_5mod = clan_membership_5
clan_membership_6mod = gtools::smartbind(data.frame(clan_membership_5mod),data.frame(FatWonderDF[which(!(FatWonderDF$profile_accession %in% clan_membership_6mod$profile_accession)),]))
clan_membership_6mod[which.duplicated(clan_membership_6mod$profile_accession),]

# clan_membership_7mod <- edit(clan_membership_8mod)
# clan_membership_8mod = distinct(clan_membership_7mod[])
# # 
# CM8t = rbind(clan_membership_6mod,filter(clan_membership, !(profile_accession %in% clan_membership_8mod$profile_accession)))
#             
#             
# CM8 = merge(CM8t,`dropcols<-`(WonderDF,c("New_Name","New_Name.x","New_Name.y","ClanID","Comp_New_Name")),by="profile_accession",all.x=T,all.y=F)
# 
# CM8$Ver = "29.09.2021"
# CM9 = rbind.fill(data.frame(CM8),data.frame(filter(`dropcols<-`(TempWonderDF,c("New_Name.x","New_Name.y","ClanID","Comp_New_Name")), !(profile_accession %in% CM8$profile_accession)),fill =T))
CM9 = clan_membership_6mod
CM9$PolyPort_Problematic = F
CM9$PolyPort_Problematic[which(CM9$profile_accession %in% PolyProts)] = T
CM9$freq = NULL
CM9$Classified = -1
CM9$Classified[which.Not.na(CM9$New_Name)] = 3
# xCM9$Classified[which.Not.na(CM9$Comment)] = 1
CM9$Classified[grep(pattern = "_fragment|CTD|NTD",CM9$New_Name)] = 1
CM11 = distinct(CM9[,c("profile_accession", "ClanID", "New_Name", "Comment","Classified","SCOPe_Function_ID","SCOPe_Class","SCOPe_PDBID","SCOPe_Desc","CDD_name","CDD_desc","CATH_name","Ver","name1","PolyPort_Problematic")])
# CM11$Fragmented[which.Not.na(CM11$New_Name)] = F
# # View(CM11[grep(pattern = "_fragment|CTD|NTD|Nterm",CM11$New_Name),])
# 
fwrite(sep = '\t',file = "CM11.tsv",x = CM11,quote = T)
# xlsx::write.xlsx(CM11,"CM11.xlsx")#


# CM10 = rbind.fill(data.frame(CM9),data.frame(filter(FatWonderDF, !(profile_accession %in% CM8$profile_accession))))

CM11 = distinct(CM9[,c("profile_accession", "ClanID", "New_Name", "Comment","SCOPe_Function_ID","SCOPe_Class","SCOPe_PDBID","SCOPe_Desc","CDD_name","CDD_desc","CATH_name","Ver","name1")])
# CM12 = distinct(edit(CM11))
# CM12 = distinct((CM12))
# CM12$New_Name[c(108722,100132,63441)] = "Lysis?"
# 
# CM12$Classified = F
# 
# CM12$Classified[which.Not.na(CM12$New_Name)] = T
# 
AllDF = merge(AllDF,CM12,by="profile_accession",all.x=T,all.y=F)
CM13 = merge(CM12,distinct(AllDF[,c("profile_accession","New_Name")]),by="profile_accession",all.x=T,all.y=F)
CM13$New_Name = CM13$New_Name.x
CM13$New_Name[which.na(CM13$New_Name.x)] = CM13$New_Name.y[which.na(CM13$New_Name.x)]
CM13$New_Name.x = NULL
CM13$New_Name.y = NULL
CM12 = distinct((CM13))

CM12$Classified[CM12$PolyPort_Problematic] = -1
CM12$Classified = 0
CM12$Classified[which.Not.na(CM12$New_Name)] = 2
CM12$Classified[grep(pattern = "?",x = CM12$New_Name,fixed = T)] = 1


# w1 = intersect(x = which.na(CM12$New_Name),y = which((pattern = "Lys",x = AllDF$matchID,fixed = T))
# AllDF$New_Name[w1] = "Lysis?"




# 
# # # # # # #
# which.duplicated(clan_membership_8mod$profile_accession)
# which.duplicated(clan_membership_6mod$profile_accession,returnValue = T)
# TempWonderDF2 = merge(WonderDF,clan_membership_3mod,by = "profile_accession",all.x=T,all.y=F)
# SharedCols(CM9,FatWonderDF) 
# 
# TempWonderDF[,p0("Comp_",sharedcols)] = F
# 
# setDF(TempWonderDF)
# TempWonderDF[,(sharedcols)] = NA
# 
# for (colnym in sharedcols) {
#   TempWonderDF[,p0("Comp_",colnym)] = (TempWonderDF[,p0(colnym,".x")] == TempWonderDF[,p0(colnym,".y")])
#   nass = intersect(which.na(TempWonderDF[,p0(colnym,".x")]),which.na(TempWonderDF[,p0(colnym,".y")]))
#   TempWonderDF[nass, p0("Comp_",colnym)] = T 
#   nass = which.na(TempWonderDF[,p0("Comp_",colnym)])
#   if(length(nass)!=0){
#     TempWonderDF[nass, p0("Comp_",colnym)] = F 
#   }
#   
#   TempWonderDF[,colnym] = TempWonderDF[,p0(colnym,".x")]
#   w1 = which.na(TempWonderDF[,p0(colnym,".x")])
#   w2 = which.Not.na(TempWonderDF[,p0(colnym,".y")])
#   w3 = intersect(w2,w1)
#   if(length(w3)!=0){
#     TempWonderDF[w3,colnym] = TempWonderDF[w3,p0(colnym,".y")]
#   }
# }
# 
# View(TempWonderDF[,c(sharedcols,"New_Name.x","New_Name.y")])
# 
# 
# 
# clan_membership5 = TempWonderDF %>% group_by(ClanID) %>% summarise(
#   MemsList= list(unlist(list(unique(profile_accession)))),
#   MemsChar= toString(unlist(MemsList)),
#   Nmems = lengths(MemsList),
#   NameList= list(unlist(list(unique(na.omit(New_Name))))),
#   NamesChar= toString(unlist(NameList)),
#   Nnames = lengths(NameList)#,
#   # New_Name = New_Name
# )
# 
# clan_membership6 = `dropcols<-`(clan_membership5,c("MemsList","NameList"))
# clan_membership6 = TabUnroll(dfdt = clan_membership6,sep = ", ",colnym = "MemsChar",NewColNym = "profile_accession")
# clan_membership6 = merge(clan_membership6,TempWonderDF[,c("New_Name","profile_accession")],by = "profile_accession", all.x=T,all.y=F)
# clan_membership6 = clan_membership6[-which.na(clan_membership6$ClanID),]
# xlsx::write.xlsx(clan_membership6,"clan_membership6.xlsx")# 
# xlsx::write.xlsx(TempWonderDF,"TempWonderDF.xlsx")# 
# 
# 
# w1 = which.na(Profs27092$New_Name)
# w2 = intersect(Profs27092$profile_accession[w1], TempWonderDF$profile_accession[which.Not.na(TempWonderDF$New_Name)])
# sum(as.integer(na.omit(HMMstats15$freq[which(HMMstats15$profile_accession %in% w2)])))



# {
# TempProfs = fread("temp_ProfsFrom_ProfCounts0209tax_20210917.tsv",colClasses = "character",header = T,fill=T,na.strings = "" )
# View(TempProfs[which.duplicated(TempProfs$profile_accession),])
# # TempProfs=merge(TempProfs,tmmp[w1,],by="profile_name",all.x=T,all.y=F) 
# # TempProfs = edit(TempProfs)
# # tmmp2 = merge(tmmp,TempProfs,by="profile_accession",all.x=T,all.y=T)
# tmmp2 = merge(tmmp,TempProfs,by="profile_accession",all.x=T,all.y=F)
# View(tmmp2[which.duplicated(tmmp2$profile_accession),])
# View(tmmp2[which(tmmp2$New_Name.x != tmmp2$New_Name.y),])
# 
# 
#   tmmp2[,setdiff("profile_accession",SharedCols(tmmp,TempProfs))] = NA
#   compdf5 = tmmp2
#   compdf[,c("Comp_New_Name","Comp_Comment")] = F
#   setDF(compdf)
#   for (colnym in   c("New_Name","Comment")) {
#     
#     compdf5[,p0("Comp_",colnym)] = (compdf5[,p0(colnym,".x")] == compdf5[,p0(colnym,".y")])
#     nass = intersect(which.na(compdf5[,p0(colnym,".x")]),which.na(compdf5[,p0(colnym,".y")]))
#     compdf5[nass, p0("Comp_",colnym)] = T 
#     nass = which.na(compdf5[,p0("Comp_",colnym)])
#     if(length(nass)!=0){
#       compdf5[nass, p0("Comp_",colnym)] = F 
#     }
#     
#     compdf5[,colnym] = compdf5[,p0(colnym,".x")]
#     w1 = which.na(compdf5[,p0(colnym,".x")])
#     w2 = which.Not.na(compdf5[,p0(colnym,".y")])
#     w3 = intersect(w2,w1)
#     if(length(w3)!=0){
#       compdf5[w3,colnym] = compdf5[w3,p0(colnym,".y")]
#     }
#   }
#   View(compdf5[,c("profile_accession","New_Name","Comment","New_Name.x","Comment.x","New_Name.y","Comment.y")])
#   TempWonderDF
#   intersect(AllDF$profile_accession[which.na(AllDF$Analysis_Type)],TempProfs$profile_accession)
#   setdiff(colnames(TempWonderDF),SharedCols(TempWonderDF,AllDF))
#   View(CM12[,setdiff(colnames(CM12),SharedCols(CM12,AllDF))])
# }

# # # # # # #
AllDF = readRDS("../Annotations/AllDF.05102021.RDS")

IDFT$Contig[which(IDFT$NC90=="NC90_091870")]
writeXStringSet(filepath = "NC90_091870.faa",AllORFs[AllORFsInfo$ORFID[which(AllORFsInfo$seqid %in% IDFT$Contig[which(IDFT$NC90=="NC90_091870")])]])
AllORFsInfo$ORFID[which(AllORFsInfo$seqid %in% IDFT$Contig[which(IDFT$NC90=="NC90_091870")])]

AllmsDFxIDFT = merge(data.frame(IDFT),data.frame(Rename1Col(AllDF,"ND","Contig")),by="Contig",all.x=T,all.y=F)
AllmsDFxIDFT = `dropcols<-`(AllmsDFxIDFT,c("pCoverage", "evalue","end","strand","RdRp_ID","UVA","nseq","M", "eff_nseq", "relent", "info", "p_relE", "compKL"))
countDF = plyr::count(AllmsDFxIDFT[,c("Order","profile","New_Name")])
# WriteWolfTbl(AllmsDFxIDFT,"AllmsDFxIDFT.tsv")
# WriteWolfTbl(Rename1Col(AllmsDF,"seqnames","ND"),"AllmsDF.tsv")

TMP4 = plyr::count(AllmsDFxIDFT[,c("Phylum","profile")])
TMP5 = plyr::count(AllmsDFxIDFT[,c("Phylum","profile_accession")])


# # # # # # #
AllDF1 <- rbind.fill(AllDF, filter(p0002MA[, SharedCols(p0002MA,AllDF)],!(ORFID %in% AllDF$ORFID)))

fwrite0(x = AllDF1,file = "../Vfin/ORFs/AllDF.15112021.tsv",sep = '\t',quote = T)
fwrite0(x = AllDF1,file = "../Annotations/AllDF.15112021.tsv",sep = '\t',quote = T)
w1 <- wna(AllDF1$ND)
AllDF1$ND[w1] <- str_split_fixed(string = AllDF1$ORFID[w1], pattern = "__", n = 2)[,1]


