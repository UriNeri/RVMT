# Presence Absence (generic).
library(Biostrings)
library(readr)
library(data.table)
functions_file="/urigo/urineri/scripts/RNA_Viruses/v2/R_scripts/fun_ctions.r"
source(functions_file)
{ # Args. # legacy - moved to fixed params by iter number from reading a table (Presence_Absence_Params.tsv)
  args <- commandArgs(trailingOnly= T)
  Algo=as.character((args[1]))  # BlastN, DiamondX, or DC-MegaBlast.
  iter=as.numeric((args[2])) # Iteration
  pid=as.character((args[3])) #Precent iden.
  evale=as.character((args[4])) #E-value 
  scov=as.character((args[5]))  # Subject coverage. Set to a negative value if you didn't set it for the Algo run (which opts to DEF).
  qcov=as.character((args[6]))  # Query coverage. Set to a negative value if you didn't set it for the Algo run (which opts to DEF).
  PWD=as.character((args[7]))  # Current working directory/outdir. Something like "/scratch200/urineri/V2/output/" Keep the leading and trailing "/" !!!
  min_len=as.numeric((args[8])) # Minimal alignment length.
  DB=as.character((args[9])) # Subject database.
  min_bitscore=85 # only for DC-MegaBlast.
} # Args.
{ # Args for interactive jobs.
  Algo="BlastN"  # BlastN, DiamondX, or DC-MegaBlast.
  iter=9 # Iteration
  pid=75 #Precent iden. 
  evale="0.0001" #E-value
  scov=-1 # Subject coverage. Set to a negative value if you didn't set it for the Algo run (which opts to DEF).
  qcov=50  # Query coverage. Set to a negative value if you didn't set it for the Algo run (which opts to DEF).
  # PWD="/scratch200/urineri/V2/absenties/"   # Current working directory/outdir. Something like "/scratch200/urineri/V2/output/" Keep the leading and trailing "/" !!!
  min_len=250 # Minimal alignment length.
  DB="DNAvir_DB" # Subject database.
  min_bitscore=82 # only for DC-MegaBlast.
  # PWD="/media/uri/HDD1/uri/V1_hard_mode/"
} # Args for interactive jobs.
search_res_dir="/media/uri/HDD1/uri/V1_hard_mode/V2/Presence_Absence/search_results/V2_scaffolds_vs_DNAvir_DB/"
work_dir="/media/uri/HDD1/uri/V1_hard_mode/V2/V2.03/"
setwd(work_dir)
scaffolds_fasta=readDNAStringSet("../Presence_Absence/absenties/V2_scaffolds_P0.fasta")
if(Algo =="BlastN"){
colnyms=c("qseqid","sseqid","pident","length","mismatch","gapopen","qstart","qend","sstart","send","evalue","bitscore")
output_file=paste0(search_res_dir,Algo,"_V2_scaffolds_vs_",DB,"_Pid",pid,"_qcov",qcov,"_eval",evale,".tsv")
Algo_out <- try(read.delim(file =output_file, header = F, sep = "\t", col.names=colnyms,stringsAsFactors = F))
TBDiscard=Algo_out$qseqid[Algo_out$length>min_len]
Discarded_scafs=as.character(unique(TBDiscard))
absent_scafs=scaffolds_fasta[setdiff(names(scaffolds_fasta),Discarded_scafs)]
writeXStringSet(absent_scafs, paste0(search_res_dir,"V2_scaffolds_P2.fasta"))
print(paste0("For iter:   ", iter,"  the pass/total ratio is:   ",(length(absent_scafs)/(length(scaffolds_fasta))),", (passed =  ",length(absent_scafs),")"))
}



timestamp()

{  # Get matches from previous iterations.
  scaf_data=data.frame("qseqid"=names(absent_scafs),"Length"=absent_scafs@ranges@width,stringsAsFactors = F)
  scaf_data_1=merge(scaf_data,Algo_out,by=1,all.x=T,all.y=F)
  Dcolnyms=c("qseqid", "sseqid", "staxids", "stitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue","bitscore","salltitles")
  Bcolnyms=c("qseqid", "sseqid", "sscinames", "pident", "qlen", "length", "mismatch", "evalue", "bitscore", "stitle", "sskingdoms", "staxids")
  Diamond_nr_cated=paste0(PWD,"cated_diamonds.tsv")
  BlastN_nt_cated=paste0(PWD,"cated_blastn.tsv")
  Diamond_nr_out <- try(read.delim(file=Diamond_nr_cated, header = F, sep = "\t", col.names=Dcolnyms,stringsAsFactors = F))
  for (i in 1:nrow(Diamond_nr_out)){ # Heavy, but need make sure multiple hits won't get into NAs.
    Diamond_nr_out$staxids[i]=as.numeric(unlist(strsplit(Diamond_nr_out$staxids[i],split=";", fixed=T))[1])
  }
  BlastN_nt_out <- try(read.delim(file=BlastN_nt_cated, header = F, sep = "\t", col.names=Bcolnyms,stringsAsFactors = F))
  lineages_kingdonly=lineages_2019_02_20[,c(1,2)]
  colnames(lineages_kingdonly)=c("staxids","superkingdom")
  Diamond_nr_out=merge(Diamond_nr_out,lineages_kingdonly,by=c("staxids"),all.x=T,all.y=F)
  # rm(merged_Diamonds)
  scaf_data=data.frame("old_ID"=names(absent_scafs),"Length"=absent_scafs@ranges@width,stringsAsFactors = F)
  scaf_data_nr=merge(scaf_data,Diamond_nr_out,by=1,all.x=T,all.y=F)
  scaf_data_nr_no_viruses=scaf_data_nr[-which(scaf_data_nr$ID %in% scaf_data_nr$ID[which(scaf_data_nr$s)] )]
  scaf_data_nt=merge(scaf_data,BlastN_nt_out,by=1,all.x=T,all.y=F)
  # Largest_scaf = absent_scafs[which(absent_scafs@ranges@width==max(absent_scafs@ranges@width))]
  # writeXStringSet(Largest_scaf, paste0(PWD,"Largest_scaf",(iter+1),".fasta"))
  # writeXStringSet(absent_scafs["3300012699_Ga0157593_1042799"] , paste0(PWD,"3300012699_Ga0157593_1042799.fasta"))
  scaffolds_IDs_MTs_dffed <- read_delim("rmdup1200_scaffolds_IDs_MTs_dffed.txt","\t", escape_double = F, col_names = F)
  scaffolds_IDs_MTs_dffed=scaffolds_IDs_MTs_dffed[,-3]
  colnames(scaffolds_IDs_MTs_dffed)=c("new_ID","old_ID")
  scaf_data=merge(scaf_data,scaffolds_IDs_MTs_dffed,by=c("old_ID"),all.x=T,all.y=F)
  names(absent_scafs)=scaf_data$new_ID
  writeXStringSet(absent_scafs, paste0(PWD,"V2_Scaffolds.fasta"))
}                               
