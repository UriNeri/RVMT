# Information  --------------------------------------------------------------------
## Script name: KvRev1.R
## Purpose of script: Calculate some misc. numerics requested by reviewers.
## Email: uri.neri@gmail.com
# Body  --------------------------------------------------------------------

#### Regarding the identification of plant or animal pathogens - Movement protein analysis:
Novel_MP_DFT <- filter(SimonDFT, Host_Evidence == "Host related domain - Movement protein")
Novel_MP_DFT <- filter(Novel_MP_DFT, RvANI90 %ni% unlist(unname(SimonDFT[Novel == F,"RvANI90"])))
unique(plyr::count(distinct(Novel_MP_DFT[,c("Family","RvANI90")]),vars = "Family"))
len(unique(Novel_MP_DFT$RvANI90))

#### Regarding the identification of plant or animal pathogens - same RvANI90 cluster as pathogens
PathoRvANI90 <- SimonDFT[Host_Evidence == "virushostdb",c("RvANI90","Host")]
PathoRvANI90 <- PathoRvANI90$RvANI90[grep("Eukaryota",PathoRvANI90$Host)]

Novel_Pathogen_DFT <- filter(filter(SimonDFT, Novel == T), RvANI90  %in% PathoRvANI90)
len(unique(Novel_Pathogen_DFT$RvANI90))
unique(plyr::count(distinct(Novel_Pathogen_DFT[,c("Family","RvANI90")]),vars = "RvANI90"))

unique(unname(unlist(plyr::count(distinct(Novel_Pathogen_DFT[,c("Family","RvANI90")]),vars = "Family")[1])))

intersect(Novel_Pathogen_DFT$RvANI90,Novel_MP_DFT$RvANI90) #unname(unlist(plyr::count(distinct(Novel_MP_DFT[,c("Family","RvANI90")]),vars = "RvANI90")[1])),unname(unlist(plyr::count(distinct(Novel_MP_DFT[,c("Family","RvANI90")]),vars = "RvANI90")[1])))


#### Regarding False positive rate of the Contig_Set_enrichment procedure - benchmarking against a Rumen metagenome: identification of plant or animal pathogens - Movement protein analysis:
Raw_MMseqsout = fread(file = "../Revision/IDFT_vsRumenDB.tsv",header = F,col.names = c("query_name", "subject_name", "evalue", "gapopen", "pident", "nident", "q1", "q2", "qL", "p1", "p2", "pL", "ali_len", "raw", "score", "mismatch", "qcov", "tcov"),sep = "\t",nThread = THREADS)
# Oldcated <- fread(file = "../../V4/cated.tsv",header = F,col.names = c("query_name", "subject_name", "evalue", "gapopen", "pident", "nident", "q1", "q2", "qL", "p1", "p2", "pL", "ali_len", "raw", "score", "qframe", "mismatch", "qcov", "tcov"),sep = "\t",nThread = THREADS)
# Raw_MMseqsout <- rbindlist(list(Raw_MMseqsout,Oldcated),use.names = T,fill = T)
# Raw_MMseqsout$qframe <- NULL


Raw_MMseqsout = filter(Raw_MMseqsout,evalue < 0.000000001)  # Estimated to avoid capture of sequences who may not be represented in the DNA filtration pipeline.
Raw_MMseqsout = filter(Raw_MMseqsout,pident >= 95)

### Zero size table, no entries passing the criteria).
#
# Raw_MMseqsout = filter(Raw_MMseqsout,qL >=  pL) # Filter to sequences contained within ones that underwent the main DNA filtration pipeline.
# Raw_MMseqsout = filter(Raw_MMseqsout,tcov >=  0.95)  # Filtering to sequences *mostly* contained (tcov) within the HQ query contigs.
#
# Raw_MMseqsout$Source <- str_split_fixed(Raw_MMseqsout$subject_name,pattern = "_",n=2)[,1]
# Raw_MMseqsout <-filter(Raw_MMseqsout,query_name %in% SimonDFT$ND)
#
# Raw_MMseqsout$TotalMatches <- Raw_MMseqsout$pident * Raw_MMseqsout$ali_len
# setorderv(setDT(Raw_MMseqsout),c("subject_name","TotalMatches","score","evalue","ali_len"),order = c(-1,-1,-1,1,-1))
# Raw_MMseqsout <- unique(Raw_MMseqsout, incomparables=FALSE, fromLast=FALSE,by="subject_name")
# fwrite(Raw_MMseqsout,"MMseqsout_ContigsDF.tsv",sep = "\t",na = "NA",nThread = 11)
