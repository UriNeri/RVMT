# Information  --------------------------------------------------------------------
## Script name: KvHOST2s.R
## Purpose of script: Consolidating host predictions/information.
## Email: uri.neri@gmail.com
# Body  --------------------------------------------------------------------
source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
RVMT="/media/HDD1/uri/RNA_Vir_MTs/RVMT/"
RNA_Vir_MTs_dir="/media/HDD1/uri/RNA_Vir_MTs/V3/"
versiza="Wolf"
current_path=p0(RNA_Vir_MTs_dir,versiza,"/")
setwd(current_path)

THREADS=11
Memory=11800

rdrp_faa = readAAStringSet("Degapped_merged.210416f.faa")
names(rdrp_faa) = trimws(names(rdrp_faa))
rdrp_faa = rdrp_faa[intersect(names(rdrp_faa), IDFT$RdRp_ID)]

##### Get information #####
MTstats=readRDS("../Vfin/RDSs/MTstats.RDS")
IDFT.fasta=readRDS("../Vfin/RDSs/ALL_nuc_3007.RDS")
IDFT=readRDS(file = "../Vfin/RDSs/IDFT.RDS")
saveRDS(IDFT,file = "../Vfin/RDSs/IDFT.RDS")
WriteWolfTbl(IDFT, "../Vfin/IDFT.082021.tsv")
AllGFF=readRDS("../Vfin/ORFs/AllGFF.RDS") #The gff used for annotations, not for exotic genetic code assignment.
MegaDF=readRDS("../Vfin/ORFs/ORFsDF.RDS")
MiscData=readRDS("../Vfin/RDSs/MiscData.RDS")
Oddities=readRDS("../Vfin/RDSs/Oddities.RDS")
Gcodetable=fread("/media/HDD1/uri/RNA_Vir_MTs/V3/Vfin/Metadata/codetable.tsv",header = T,sep = "\t")

##### Start consolidating (Attempt #1.)#####
# IDFT$Host = "UA"

##### Based on isolate information (From Mbio18) #####
# Mbio18Hosts = fread0("../Simon/Host_prediction/Mbio18Hosts.tsv")
# Mbio18Hosts$UA = F
# Mbio18Hosts$UA[grep(pattern = "Xinzhou|Hubei|Zhejiang|Partitivirus_like|picobi|Picobi",x = Mbio18Hosts$nam,fixed = F,ignore.case = T)] = T
# Mbio18Hosts = `droprows<-`(Mbio18Hosts,which(Mbio18Hosts$host==""))
# Mbio18Hosts = `droprows<-`(Mbio18Hosts,which(Mbio18Hosts$UA))
# 
# for (i in 1:nrow(Mbio18Hosts)){
#   ix = grep(pattern = Mbio18Hosts$GenBank[i],x = IDFT$full_name,fixed = T,ignore.case = T)
#   if(length(ix)==0){next}
#   IDFT$Host[ix] = Mbio18Hosts$host[i]
#   print(p0("found host for ", IDFT$ND[ix]))
# }

##### By Non-Standard Genetic code / MMTESP #####
# IDFT4 = IDFT

MMTESPs = MTstats$MTstats$`IMG Genome ID`[grep(x = MTstats$MTstats$`Sample Name`,pattern = "MMETSP",fixed = T)]
MMETSP_Info = fread("/media/HDD1/uri/RNA_Vir_MTs/V3/Simon/MMETSP_Info.tsv",na.strings="NA")
MMETSP_Info$Lineage[wh(MMETSP_Info$IMG_UID ==3300017106 )] = "Eukaryota; Viridiplantae; Chlorophyta; core chlorophytes; Chlorodendrophyceae; Chlorodendrales; Chlorodendraceae"
MMETSP_Info$GeneticCode.GCName[wh(MMETSP_Info$IMG_UID ==3300017106 )] = "Standard"

MMETSP_Info = distinct(MMETSP_Info[,c("Lineage","IMG_UID","GeneticCode.GCName")])
colnames(MMETSP_Info) = c("Host","IMG.Genome.ID","Host_GeneticCode_Name")
MMETSP_Info$IMG.Genome.ID = as.char(MMETSP_Info$IMG.Genome.ID)
setdiff(MMTESPs,MMETSP_Info$IMG.Genome.ID)

MMETSP_Info = filter(MMETSP_Info,IMG.Genome.ID %in% na.omit(IDFT3$Source))
MMETSP_Info$Host = gsub(pattern =  "cellular organisms; " ,replacement = "",x = MMETSP_Info$Host,fixed = T)
for (i in 1:nrow(MMETSP_Info)){
  ix = wh(IDFT3$Source == MMETSP_Info$IMG.Genome.ID[i])
  if(length(ix)==0){next}
  # IDFT4$Host[ix] = MMETSP_Info$Host[i]
  print(p0("found host for ", IDFT3$ND[ix]))
  print(MMETSP_Info$Host_GeneticCode_Name[i])
  if(MMETSP_Info$Host_GeneticCode_Name[i] =="Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear"){
    IDFT3$Host_Evidence[ix] <- "MMETSP Host genetic code matches the virus code"
    IDFT3$Host[ix] <- MMETSP_Info$Host[i]
  }
}
IDFT = IDFT3
IDFT$node.labe=NULL
IDFT <- merge(setDF(IDFT),tmpdf,by="node",all.x=T,all.y=F)
dupIDtts <- IDFT[whd(IDFT$ND),]
IDFT <- IDFT[-whd(IDFT$ND),]
IDFT <- distinct(IDFT)
# 
fwrite(IDFT,"../Wolf/IDDT8.15112021.tsv",na = "",sep = '\t')

##### 11.11.2021 #####
IDFT$Host_Evidence <- NA
IDFT$Host <- NA
IDFT4 <- IDFT
# IDFT4[,c("node","node.label")] = NA
# IDFT4[grep(pattern = "Lvl [0123]",x = IDFT4$AfLvl),c("node","node.label")] <- IDFT[grep(pattern = "Lvl [0123]",x = IDFT4$AfLvl),c("node","node.label")] 

##### Based on VirusHostDB (fetched 08.08.2021) #####
VHDB = fread0("/media/HDD1/uri/DBs/VirHostDB/virushostdb.daily.Riboviria.tsv")

VHDB$UA = F
VHDB$UA[grep(pattern = "ShiM-2016|Partitivirus_like|picobi|Picobi",x = VHDB$`virus lineage`,fixed = F,ignore.case = T)] = T
VHDB$UA[grep(pattern = "ShiM|Partitivirus_like|picobi|Picobi",x = VHDB$`virus name`,fixed = T,ignore.case = T)] = T

VHDB = `droprows<-`(VHDB,which(VHDB$UA))
VHDB = `droprows<-`(VHDB,which(VHDB$`host lineage`==""))

VHDB$refIDs = vector('list',nrow(VHDB))   # List place holder 1
VHDB$refIDs = apply(X =VHDB,MARGIN = 1,FUN = function(x) str_split(string = x["refseq id"],pattern = ", ",n = Inf,simplify = F) )

MiniVHDB <- distinct(VHDB[,c("refseq id","refIDs","host lineage")])

for (i in 1:nrow(MiniVHDB)){
  iy = unlist(MiniVHDB$refIDs[i])
  ix = NA
  for (iz in 1:length(iy)){
    ix = c(ix,grep(pattern = iy[iz],x = IDFT4$Full_name,fixed = T,ignore.case = T))
  }
  ix = ix[wana(ix)]
  if(length(ix)==0){next}
  IDFT4$Host[ix] = MiniVHDB$`host lineage`[i]
  IDFT4$Host_Evidence[ix] <- "virushostdb"
  
  print(p0("found host for ", IDFT4$ND[ix]))
}

##### By affiliation to known-host clades #####
w1 <- intersect(wna(IDFT4$Host),c(wh(IDFT4$Class %in% c("Leviviricetes")),wh(IDFT4$Family %in% c("Cystoviridae"))))
IDFT4$Host[w1] = "Bacteria"
IDFT4$Host_Evidence[w1] <- "Host known clade - Leviviricetes OR Cystoviridae"

# # w1 = which(IDFT4$Family %in% c("Picobirnaviridae"))
# # IDFT4$Host[w1] = "Host Disputed"

# Movement protein info:
w1 = grep(pattern = "MP_",x = AllDF$New_Name,fixed=T)
Putative_Movement_Encoding = unique(na.omit(AllDF$ND[w1]))
Putative_Movement_Encoding_Nodes = unique(na.omit(IDFT4$node[which(IDFT4$ND %in% Putative_Movement_Encoding)]))
w1 <- intersect(wna(IDFT4$Host),wh(IDFT4$node %in% Putative_Movement_Encoding_Nodes))
IDFT4$Host[w1] <- "Eukaryota"
IDFT4$Host_Evidence[w1] <- "Host related domain - Movement protein"


##### By presence of Bacteriophage like/related lysis gene #####
Lysers = AllDF$ND[grep(pattern = "Lysis",x = AllDF$New_Name,ignore.case = T)]
Lysers_nodes <- unique(na.omit(IDFT4$node[which(IDFT4$ND %in% Lysers)]))
w1 <- intersect(wna(IDFT4$Host),wh(IDFT4$node %in% Lysers_nodes))
w2 <- intersect(w1, c(
  wh(IDFT4$Class=="Duplopiviricetes"),
  wh(IDFT4$Phylum %in%c("p.0002.base-Kitrino", "p.0001.base-Kitrino")),
  wna(IDFT4$Phylum)))

IDFT4$Host[w2] = "Bacteria"
IDFT4$Host_Evidence[w2] <- "Putative Host related domain - bacteriolytic protein"

##### By CRISPR spacer matches #####
CRISPR = wh(IDFT4$node %in% unique(IDFT4$node[c(wana(IDFT4$Hit.s.),wana(IDFT4$Type.hit))]))

w2 <- intersect(wna(IDFT4$Host),CRISPR)
w3 <- intersect(w1, c(
  wh(IDFT4$Class=="Duplopiviricetes"),
  wh(IDFT4$Phylum %in%c("p.0002.base-Kitrino", "p.0001.base-Kitrino")),
  wna(IDFT4$Phylum)))

IDFT4$Host[w2] = "Bacteria"
IDFT4$Host_Evidence[w2] <- "Putative CRISPR spacer match"

##### Literature (manual curation by Valerian) #####
IDFT4$Host[wh(IDFT4$RID == "rv20_v2_2738")] = "Eukaryota; Viridiplantae; Chlorophyta; Ulvophyceae; TCBD clade; Bryopsidales; Bryopsidineae; Bryopsidaceae; Bryopsis; Chloroplast"
IDFT4$Host_Evidence[wh(IDFT4$RID == "rv20_v2_2738")] = "PMID:12777056"
IDFT4$Host[wh(IDFT4$RID %in% c("rv20_v2_2876","rv20_v2_4192"))] = "Eukaryota; Viridiplantae; Chlorophyta; Ulvophyceae; TCBD clade; Bryopsidales; Bryopsidineae; Bryopsidaceae; Bryopsis; Mitochondria"
IDFT4$Host_Evidence[wh(IDFT4$RID %in% c("rv20_v2_2876","rv20_v2_4192"))] = "PMID:9526504, PMID:24999047"

#### Intermission ####
IDFT4
# IDFT4 <- IDFT
IDFT4[,c("node","node.label")] = NA
IDFT4[,c("node","node.label")] <- IDFT[,c("node","node.label")] 
IDFT4[grep(pattern = "Lvl [0123]",x = IDFT4$AfLvl),c("node","node.label")] <- IDFT[grep(pattern = "Lvl [0123]",x = IDFT4$AfLvl),c("node","node.label")] 

IDFT4[,c("Phylum","Class","Order","Family","Genus")] = NA
IDFT4[grep(pattern = "Lvl [0123]",x = IDFT4$AfLvl),c("Phylum","Class","Order","Family","Genus")] = 
  IDFT[grep(pattern = "Lvl [0123]",x = IDFT4$AfLvl),c("Phylum","Class","Order","Family","Genus")]   
  
fwrite(IDFT4)




# unique(IDFT$Hit.s.[wana(IDFT$Hit.s.)])
# 
# 
# New_CRISPR=fread("/media/HDD1/uri/RNA_Vir_MTs/V3/Simon/New_CRISPR.tsv",header = T,sep = "\t")
# # setdiff(New_CRISPR$full_name,IDFT$full_name)
# # IDFT4 = merge(IDFT,New_CRISPR,by="full_name",all.x=T,all.y=F)
# # IDFT = IDFT4
# # tmpcount = plyr::count(IDFT[,c("Type hit" ,"Family")])
# # IDFT$Note = NA
# IDFT$Host[which(IDFT$NC90 == "c210416.00399")] = "Bacteria; Chloroflexi; Chloroflexia; Chloroflexales; Roseiflexaceae; Roseiflexus"
# IDFT$Note[which(IDFT$NC90 == "c210416.00399")] = "Host assignment via CRISPR match."

##### By 90% ID cluster membership #####
C9H = distinct(IDFT[,c("NC90","Host")])

C9H = `droprows<-`(C9H,which(C9H$Host == "UA"))
C9H$Host[which.duplicated(C9H$NC90,returnalldups = T,returnValue = F)] = "Eukaryota"

C9H = distinct(C9H)
IDFT4 = IDFT
IDFT4$Host = NULL
IDFT4 = merge(IDFT4,C9H, by="NC90",all.x=T,all.y=F)
IDFT4$Host[which.na(IDFT4$Host)] = "UA"
IDFT = IDFT4


##### Misc ####
# IDFT$source[grep(pattern = "lcl|",x = IDFT$source,fixed = T)]="RefSeq/GenBank"
# IDFT$Host="UA"
# IDFT$Host[which(IDFT$Class=="Leviviricetes")] = "Bacteria"
# IDFT$Host[which(IDFT$Order=="Mindivirales")]="Bacteria"
# IDFT$Host[which(IDFT$Family=="Picobirnaviridae")]="Host_Disputed"
# IDFT$Host[which(IDFT$Family=="Partitiviridae")]="Host_Disputed"
# IDFT$Host[which(IDFT$Host=="UA")]="Euk" #### Assuming assumptions !!!
# plyr::count(IDFT$Host)
# 
# 
# 
# 
# ##### 06.12.2021 - RBS (yet again) ##### 
# tmpIDFT <- IDFT
# tmpIDFT <- tmpIDFT[wana(tmpIDFT$RBS),]
# tmpIDFT <- tmpIDFT[wh(tmpIDFT$RBS != "?"),]
# 
# Clades2inspect <- c("Picobirnaviridae","Leviviricetes","Tymovirales","Ghabrivirales","Chuviridae","Cystoviridae",)
# 
# 
# 
# 
# 
# 
# 



