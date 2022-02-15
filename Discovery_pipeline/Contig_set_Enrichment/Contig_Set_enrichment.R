# Information  --------------------------------------------------------------------
## Script name: Contig_Set_enrichment.R
## Email: uri.neri@gmail.com
# Body  --------------------------------------------------------------------
source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
RVMT <- "/media/HDD1/uri/RNA_Vir_MTs/RVMT/"
RNA_Vir_MTs_dir <- "/media/HDD1/uri/RNA_Vir_MTs/V3/"
versiza <- "Simon"
current_path <- p0(RNA_Vir_MTs_dir, versiza, "/")
setwd(current_path)
library(ape)
library(castor)
library(phangorn) 
library(Cairo)

THREADS <- 11
Memory <- 11800
Raw_MMseqsout = fread(file = "SimonDFT_vs_MTs.tsv",header = F,col.names = c("query_name", "subject_name", "evalue", "gapopen", "pident", "nident", "q1", "q2", "qL", "p1", "p2", "pL", "ali_len", "raw", "score", "mismatch", "qcov", "tcov"),sep = "\t",nThread = THREADS)
# Oldcated <- fread(file = "../../V4/cated.tsv",header = F,col.names = c("query_name", "subject_name", "evalue", "gapopen", "pident", "nident", "q1", "q2", "qL", "p1", "p2", "pL", "ali_len", "raw", "score", "qframe", "mismatch", "qcov", "tcov"),sep = "\t",nThread = THREADS)
# Raw_MMseqsout <- rbindlist(list(Raw_MMseqsout,Oldcated),use.names = T,fill = T)
# Raw_MMseqsout$qframe <- NULL


Raw_MMseqsout = filter(Raw_MMseqsout,evalue < 0.000000001)  # Estimated to avoid capture of sequences who may not be represented in the DNA filtration pipeline.
Raw_MMseqsout = filter(Raw_MMseqsout,pident >= 95) 
Raw_MMseqsout = filter(Raw_MMseqsout,qL >=  pL) # Filter to sequences contained within ones that underwent the main DNA filtration pipeline.
Raw_MMseqsout = filter(Raw_MMseqsout,tcov >=  0.95)  # Filtering to sequences *mostly* contained (tcov) within the HQ query contigs.

Raw_MMseqsout$Source <- str_split_fixed(Raw_MMseqsout$subject_name,pattern = "_",n=2)[,1]
Raw_MMseqsout <-filter(Raw_MMseqsout,query_name %in% SimonDFT$ND) 

Raw_MMseqsout$TotalMatches <- Raw_MMseqsout$pident * Raw_MMseqsout$ali_len
setorderv(setDT(Raw_MMseqsout),c("subject_name","TotalMatches","score","evalue","ali_len"),order = c(-1,-1,-1,1,-1))
Raw_MMseqsout <- unique(Raw_MMseqsout, incomparables=FALSE, fromLast=FALSE,by="subject_name")
fwrite(Raw_MMseqsout,"MMseqsout_ContigsDF.tsv",sep = "\t",na = "NA",nThread = 11)
# 
# 
# 
# 
# 
# w1 <- grep(pattern = "New_nucc_id_",x =Raw_MMseqsout$query_name,fixed = T )
# Raw_MMseqsout$query_name[w1] <- gsub(pattern = "New_nucc_id_",replacement = "",x = Raw_MMseqsout$query_name[w1])
# Raw_MMseqsout$query_name[w1] <- p0("ND_",pad(str = Raw_MMseqsout$query_name[w1],pad = "0",width = 6,side = "left" ))
# 
# 
# 
# tmpdf <- distinct(Raw_MMseqsout[,c("query_name","qL")])
# tmpdf <- (tmpdf[whd(tmpdf$query_name),])
# tmpdf$Length <-IDFT.fasta[tmpdf$query_name]@ranges@width
# tmpdf <- tmpdf[wh(tmpdf$qL == tmpdf$Length),]
# tmpdf2 <- Raw_MMseqsout
# 
# for(ix in 1:nrow(tmpdf)){
#   remo2 <-  intersect(wh(tmpdf2$qL != tmpdf$Length[ix]),wh(tmpdf2$query_name == tmpdf$query_name[ix]))
#   tmpdf2 <- tmpdf2[-remo2,]
# }
# Raw_MMseqsout <- tmpdf2 
# 
# 
tmpdf4Simon <- distinct(NovoDFT[wh(NovoDFT$RID == NovoDFT$RCR90),c("RID","RCR90","ND")])
tmpdf4Simon <- tmpdf4Simon[wana(tmpdf4Simon$RID),]
WriteWolfTbl(tmpdf4Simon,"tmpdf4Simon.tsv")


SimonDFT2