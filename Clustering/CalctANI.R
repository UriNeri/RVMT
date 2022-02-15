# Information -------------------------------------------------------------------- 
# Script name: CalcANI.R 
# Purpose of script: Calculation of total average nucleic identity (modified - "Rv"ANI)
# Author: Uri Neri 
# Date Created:  08.11.2021
# Email: uri.neri@gmail.com 

# Body --------------------------------------------------------------------

# # Bash::
# module load miniconda/miniconda3-4.7.12-environmentally
# conda activate /urigo/DBs/envs/NeriEnv2
# cdtome
# R
# 
source("./basicf.R")
setwd("./AvA/tANI")

THREADS <- 46
Memory <- 11800

CallMMseqs <- function(input_path = "input.fasta", threads = 12, output_path = "AVA.mm.out", min_id = 0.75, max_evalue = 0.000001, storage_limit = "1T"){
  cmd = paste("mmseqs easy-search", input_path, input_path, output_path, " tmp -e", max_evalue, " -a true --add-self-matches true -s 3 --search-type 3 --max-seqs 5000000 --min-aln-len 100 --format-output query,target,nident,pident,qstart,qend,tstart,tend,tlen,qlen,alnlen,evalue,bits,tcov,qcov ", sep=" ")
  GenericRunner(command = cmd, param =  list("threads" = threads, "min-seq-id" = min_id, "disk-space-limit" = storage_limit) ,RunInTmp = T, FlgType = "--", AsType = " ")
}

ReadAVA <- function(File = "AVA.mm.out", threads = THREADS, ...){
  fread(input = File, quote = '"', header = F, sep = "\t", na.strings = "", 
        col.names = c("qid", "tid", "nid", "pid", "q1", "q2", "t1", "t2", "tL", "qL", "alnlen", "evalue", "bits", "tcov", "qcov"),
        # drop = c(13,14,3),
        nThread = threads,
        data.table = F,
        ...)
}


IDFT <- fread("IDFT.28",sep = "\t",header = T,quote = '"',na.strings = "NA")
IDFT.fasta <- readDNAStringSet("IDFT.fna")
contigsdf <- fread("ContigDF.27.01.2022.tsv")
AllContigs = readDNAStringSet("All.101021.fna") # 2.5Mil contig set.
AllContigs <- AllContigs[intersect(AllContigs@ranges@NAMES,contigsdf$ND)]
AllContigs <- c(AllContigs,IDFT.fasta)
AllContigs <- AllContigs[unique(names(AllContigs))]

AllvsALL <- ReadAVA("All.3101.blout")
# As the file name implies, I pre-culled weak hits using AWK to be able to work on the matrix and graph in-memory.
#  Awk's line by line approach is more similar to pythonic solutions that iterate over the blast/mmseqs output. 
#  
# w1 <- intersect(intersect(wh(AllvsALL$pid>= 75),wh(0.0001 >= AllvsALL$evalue)),wh(AllvsALL$bits >= 25))
w1 <- intersect(wh(AllvsALL$qid %in% names(AllContigs)),wh(AllvsALL$tid %in% names(AllContigs)))
len(w1) # 1575606558
nrow(AllvsALL) # 1837542584
AllvsALL <- AllvsALL[w1,]
gc()
fwrite(AllvsALL,"All.3101.blout",sep="\t",na = "NA",nThread = THREADS)

#
AllvsALL <- fread("All.3101.blout",sep="\t",na.strings = "NA",nThread = THREADS)
# gc()


# Cull each {qid,tid} pair to their best available alignment:
setDT(AllvsALL)
setorderv(AllvsALL,c("qid", "tid",  "bits"), order = c(-1,-1,-1))
AllvsALL = unique(AllvsALL, incomparables = FALSE, fromLast = FALSE, by = c("qid","tid"))
nrow(AllvsALL) # 1575767346

gc()
# Calculate average nucleic identity:
AllvsALL[, ANI := (pid * alnlen) / pmin(qL, tL)]

# Set alignment fraction by the smaller:
AllvsALL[, AF := (pmin(qcov, tcov))] # # AllVsAll[, AF := (pmin(qL, tL) / pmax(qL, tL))]

gc()

# Prune hits:
w0 <- wh(AllvsALL$ANI >= 90)
ww0 <- wh(AllvsALL$AF >= 0.9)
w1 <- intersect(w0,ww0)
rm(w0,ww0)
gc()

MicroAVA <- AllvsALL[w1,c("qid","tid")]
fwrite(MicroAVA,"AvA_ANI90_AF90.3101.AB",sep="\t",na = "NA",nThread = THREADS)

gc()

# Reinsert relatively short but perfect terminal matches (a.k.a the "Rv Correction"):
w0 <- wh(AllvsALL$pid >= 99)
ww0 <- wh(AllvsALL$alnlen >= 150)
w3 <- intersect(w0,ww0)
rm(w0,ww0)
gc()

w2 <- unique(setdiff(w3,w1))
rm(w3)
gc()

PerfectAVA <-(AllvsALL[w2,])

RevMatches <- wh(PerfectAVA$q1 > PerfectAVA$q2)

Query5Covered <- c(wh(PerfectAVA$q2 == PerfectAVA$qL),RevMatches[wh(PerfectAVA$q1[RevMatches] == PerfectAVA$qL[RevMatches])])
Query3Covered <- c(wh(PerfectAVA$q1 == 1),RevMatches[wh(PerfectAVA$q2[RevMatches] == 1)])
TerminalAlig <- unique(c(Query5Covered,Query3Covered))
rm(Query3Covered,Query5Covered)
gc()
PerfectAVA <- PerfectAVA[TerminalAlig,c("qid","tid")]
fwrite(PerfectAVA,"AvA_Perfect_terminal_alignments.AB",sep="\t",na = "NA",nThread = THREADS)

rm(AllvsALL)

# nrow(PerfectAVA)- len(TerminalAlig) 

Clinks4Gava <- rbindlist(list(MicroAVA,PerfectAVA),fill=TRUE,idcol=T,use.names = T)
Clinks4GavaPruned = Clinks4Gava[intersect(wh(Clinks4Gava$qid %in% IDDT7$ND),wh(Clinks4Gava$tid %in% IDDT7$ND)),]
gc()

# Format as an igraph object:
gAvA <- igraph::graph_from_edgelist(as.matrix(Clinks4Gava[,c("qid","tid")]), directed = F)

Cpm1 <- igraph::clusters(gAvA, mode="weak")
Cpm1$no
# summary(Cpm1$csize)
# # ggplot2: ::histogram(Cpm1$csize,)
# # geom_histogram(mapping = Cpm1$csize)
pdf("Cpm1Hist.All.TerminalAligConnected.pdf", 
    width = 10, 
    height = 10)  
ggplot(data.frame("size" = Cpm1$csize),aes(size)) +  geom_histogram(bins = 1500) + scale_y_log10()
dev.off()

gAvA2 <- igraph::graph_from_edgelist(as.matrix(AllvsALL[AF >= 0.95][ANI >= 75][bits >= 95][,c("qid","tid")]), directed = F)
Cpm2 <- igraph::clusters(gAvA2, mode="weak")

# ggplot(data.frame("size" = Cpm2$csize),aes(size)) +  geom_histogram(bins = 150) + scale_y_log10()
# summary(Cpm2$csize)


tmpdf <- data.frame("ND" = names(Cpm1$membership), ANI90 = unname(Cpm1$membership))
tmpdf$ANI90 <- p0("RvANI90_",pad(tmpdf$ANI90,pad = "0",width = 1+max(nchar(tmpdf$ANI90)),side = "left",use_length = F))
contigsdf <- merge(contigsdf, tmpdf, by="ND", all.x=T, all.y=F)
contigsdf <- merge(contigsdf, AllContigsDF[,c("ND","Length")], by="ND", all.x=T, all.y=F)

## Merge into IDFT:
IDDT8 <- merge(IDFT, tmpdf3, by="ND", all.x=T, all.y=F)
IDDT8 <- merge(IDDT8, tmpdf2, by="ND", all.x=T, all.y=F)
IDDT8 <- merge(IDDT8, tmpdf, by="ND", all.x=T, all.y=F)
saveRDS(IDDT8,"../Wolf/IDDT7.RDS")
saveRDS(IDDT8,"../Vfin/RDSs/IDDT7.RDS")
fwrite0(x = IDDT8,file = "../Wolf/IDDT7.tbl")
fwrite0(x = IDDT8,file = "../Vfin/Metadata/IDDT7.tbl")

fwrite0(x = distinct(IDDT8[,c("ND","ANI75","ANI90","ANI95")]),file = "ANI.tbl")
fwrite(compress = )
fwrite0(x = AllvsALL,file = "IDFT.AvA.tbl")


#### Obsolete: ####
# CalcAF <- function(AVA){
#   # AF = ∑(Length of the shorter fragment)/(Length of the Query Genome)
#   # Doing the sigma outside, this is only per the query side.
#   #AF <- 
#   sum(min(AVA[c("qL", "tL")]) / (AVA["qL"]))
#   # return(AF)
# }
# 
# CalcANI <- function(AVA){
#   # ANI = ∑(ID%*Length of Alignment)/(∑(Length of the shorter fragment))
#   ANI <- sum(AVA["pident"] * AVA["alnlen"]) / min(AVA[c("qL","tL")])
# }
# 
# CalcTANI <- function(AF, ANI){
#   # tANI = -ln(AF*ANI)
#   tANI <- -ln(AF * ANI)
# }


