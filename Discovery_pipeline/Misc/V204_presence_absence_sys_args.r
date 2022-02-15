#V204 - "V2-like" filtering / Presence-Absence for Sequence sets regardless if they're from studies with MGs.
##### Set env ##### 
library(Biostrings)
library(readr)
library(data.table)

##### Arguments  ##### 
{ # Args. # legacy - moved to fixed params by iter number from reading a table (Presence_Absence_Params.tsv)
  # functions_file="/urigo/urineri/scripts/RNA_Viruses/Neri_scripts/RdRp_search/Rscripts/V2_function.r"
  # source(functions_file)
  # args <- commandArgs(trailingOnly= T)
  # Algo=as.character((args[1]))  #Search tool used - MMseqs2, BlastN, DiamondX, or DC-MegaBlast.
  # pid=as.character((args[2])) #Precent iden.
  # evale=as.character((args[3])) #E-value 
  # PWD=as.character((args[4]))  # Current working directory/outdir. Something like "/scratch200/urineri/V2/output/" Keep the leading and trailing "/" !!!
  # min_len=as.numeric((args[5])) # Minimal alignment length.
  # Algo_outputP=as.character((args[6])) # Minimal alignment length.
  # Input_fasta=as.character((args[7])) # Minimal alignment length.
  search_res_dir="/scartch200/urineri/V2/MTs/step_8/"
  # s=as.character((args[8]))
  # stud=as.character((args[9]))
  # Algo_outputP=file.path(search_res_dir,s,stud,paste0("MMseqs2_JGI_set_vs_",stud,".tsv"))
  # fasta_dir="/scartch200/urineri/V2/MTs/results_dir/fasta_files/"
  # Input_fasta=file.path(fasta_dir,"Putative_RNA_Viruses.fna")
  Non_disqualifyiers_file="/urigo/urineri/RNA_Viruses/acc/ND_list.txt"
  # functions_file="/urigo/urineri/scripts/RNA_Viruses/v2/R_scripts/fun_ctions.r"
  lineages <- fread("/urigo/urineri/RNA_Viruses/acc/lineages_merged.tsv") # Get all for now, later maybe keep only Kingdome
  
} # Args.
{ # Args for interactive jobs.
  # functions_file="/media/uri/HDD1/uri/V1_hard_mode/V2/Neri_scripts/RdRp_search/Rscripts/V2_function.r"
  # source(functions_file)
  Algo="DiamondX"  # BlastN, DiamondX, or DC-MegaBlast.
  pid=55 #Precent iden. 
  evale=0.0001 #E-value
  PWD="/media/uri/HDD1/uri/V1_hard_mode/"
  # PWD="/scratch200/urineri/V2/MTs/"   # Current working directory/outdir. Something like "/scratch200/urineri/V2/output/" Keep the leading and trailing "/" !!!
  min_len=30# Minimal alignment length.
  # min_bitscore=82 # only for DC-MegaBlast.
  # search_res_dir="/media/uri/HDD1/uri/V1_hard_mode/V2/Presence_Absence/V204/search_results/"
  # stud="Lab_enriched_sorghum-adapted_microbial_communities_from_Joint_BioEnergy_Institute__Emeryville__California__United_States"
  # Algo_output_file=file.path(search_res_dir,"cated")
  # fasta_dir="/media/uri/HDD1/uri/V1_hard_mode/V2/Presence_Absence/V204/fasta_files/"
  # Input_fasta_file=file.path(fasta_dir,"JGI_putative_RNA_Viruses.fna")
  s="MGs"
  Non_disqualifyiers_file="/media/uri/HDD1/uri/V1_hard_mode/V2/acc/ND_list.txt"
  lineages <- fread(paste0(PWD,"/V2/acc/lineages_merged.tsv"),na.strings = "NA") # Get all for now, later maybe keep only Kingdome
  Manual_TN=c("3300030607_Ga0247615_10000884","3300028184_Ga0256796_1013561", "3300028179_Ga0256903_1011018","3300028171_Ga0256791_1021197","3300026493_Ga0256786_1026851")
  
}
##### Set consts.  ##### 
Cellular_nms= c("Archaea","Eukaryota","Bacteria")
Red_Herrings=c("virus","phage","capsid","coat","tail","LysM","lysis", "viral","vector","AngRem","PRSV","ds-RNA element","dsRNA element","Bromo_MP","TBSV_P22","Tombus_P19","Arena_ncap_C","Peptidase_A21","Peptidase_A6")
constructs_taxIDs=c("32630")
BLASTn_vNT_colnyms=c("qseqid", "sseqid", "sscinames", "pident", "qlen", "alnlen", "mismatch", "evalue", "bits", "stitle", "sskingdoms", "tax_id")
MMseqs_vDNA_colnyms=c("qseqid","sseqid","evalue","pident","qstart","qend","qlen","tstart","tend","tlen","alnlen","raw","bits","mismatch")
BLASTN_vIMGVR_colnyms=c("qseqid","sseqid","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")
Diamond_vNR_colnyms=c("qseqid", "sseqid", "tax_id", "stitle", "pident", "alnlen", "mismatch", "gapopen", "qstart", "qend", "tstart", "tend", "evalue","bits","salltitles")
Diamond_vMGs_colnyms=c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "evalue", "bits", "alnlen", "pident", "mismatch", "slen")

# MMseqs_vPfam_colnyms
Non_disqualifyiers=scan(Non_disqualifyiers_file, what="", sep="\n")
manual_FPs=scan("presence_absence/new_manual_fps.txt", what="", sep="\n")

# Euka_RdRps_MMseqs2=scan(paste0(PWD,"V2/Presence_Absence/V204/search_results/Euka_RdRps.txt"), what="", sep="\n")
# FPs=union(Euka_RdRps_MMseqs2,Manual_TN)
setwd(PWD)
# setwd("step_8")
setwd("V2/V2.04/")


##### Helper functions  ##### 
filter_by_LeP <- function(Algo_out=Algo_out,min_len=min_len,pid=pid,evale=evale){
  legnties=which(Algo_out$alnlen>=min_len)
  valars=which(Algo_out$evalue<=evale)
  persenti=which(Algo_out$pident>=pid)
  keepers=intersect(intersect(persenti,valars),legnties)
  Algo_out=Algo_out[keepers,]
  return(Algo_out)
}

##### input vs DNAome ##### 
# step=4
# s="NR"
# c="0"
# 
# # setwd(file.path(PWD,paste0("step_",step),paste0("C",c)))
# # Input_fasta_file=paste0("step_",(step-1),"_filtered_ml1200.",c,".fasta")
# Input_fasta_file="DNAVR_fitered.fasta"
# # Algo_output_file=paste0("MMseqs2_C",c,"_vs_",s,".tsv")
# Algo_output_file="DiamondX_vs_NR.tsv"
# Input_fasta=readDNAStringSet(Input_fasta_file)
MMseqs2_vMGs=fread(Algo_output_file,col.names = MMseqs_vDNA_colnyms,stringsAsFactors = F)
MMseqs2_vMGs=filter_by_LeP(Algo_out=MMseqs2_vMGs,min_len=min_len,pid=pid,evale=evale)
# # min_dt=data.frame(MMseqs2_vMGs[ , max(bits), by = qseqid])
# # colnames(min_dt)=c("qseqid","bits")
# # MMseqs2_vMGs=merge(min_dt,MMseqs2_vMGs,by=c("qseqid","bits"),all.x=T,all.y=F)
# # MMseqs2_vMGs=MMseqs2_vMGs[-which(duplicated(MMseqs2_vMGs$qseqid)),]
removies=unique(MMseqs2_vMGs$qseqid)
# # removies=intersect(unique(MMseqs2_vMGs$qseqid),names(Input_fasta))
keepers=setdiff(names(Input_fasta),removies)
output_fasta=Input_fasta[keepers]
writeXStringSet(output_fasta,paste0(s,"_fitered_C",c,".fasta"))

##### input vs IMG/VR (DNA only) ##### 
# step=4
# s="IMG_VR_DNA"
# c="1"
# 
# setwd(file.path(PWD,paste0("step_",step),paste0("C",c)))
# Input_fasta_file=paste0("step_",(step-1),"_filtered_ml1200.",c,".fasta")
# Algo_output_file=paste0("MMseqs2_C",c,"_vs_",s,".tsv")
# Input_fasta=readDNAStringSet(Input_fasta_file)
# 
# BLASTn_VR=fread(file =paste0(search_res_dir,"/","/BlastN_JGI_putative_RNA_Viruses_vs_DNAvir_DB_Pid75_eval0.0001.tsv"), header = F, sep = "\t", col.names=BLASTN_vIMGVR_colnyms,stringsAsFactors = F)
# legnties=which(BLASTn_VR$alnlen>=min_len)
# valars=which(BLASTn_VR$evalue<=evale)
# persenti=which(BLASTn_VR$pident>=(pid*100))
# keepers=intersect(intersect(persenti,valars),legnties)
# BLASTn_VR=BLASTn_VR[keepers,]
# FPs=unique(append(FPs,BLASTn_VR$qseqid))
# 
##### input vs MGs (DiamondX) ##### 
Input_fasta_file="FPs_filtered.fasta"
Input_fasta=readDNAStringSet(Input_fasta_file)
Algo_out=fread("cated_step_7_vs_DNAome.tsv",stringsAsFactors = F,col.names = Diamond_vMGs_colnyms)
min_len=30
pid=30
evale=0.0001
Algo_out=filter_by_LeP(Algo_out=Algo_out,min_len=min_len,pid=pid,evale=evale)
min_dt=data.frame(Algo_out[ , max(bits), by = qseqid])
colnames(min_dt)=c("qseqid","bits")
Algo_out=merge(min_dt,Algo_out,by=c("qseqid","bits"),all.x=T,all.y=F)
removies=unique(Algo_out$qseqid)
keepers=setdiff(names(Input_fasta),removies)
output_fasta=Input_fasta[keepers]
writeXStringSet(output_fasta,"DiamondX_MG_filtered.fasta")
##### input vs NCBI's NR ##### 
Algo="DiamondX"
Input_fasta_file="DNAVR_filtered.fasta"
Input_fasta=readDNAStringSet(Input_fasta_file)
# Algo_output_file=paste0("MMseqs2_C",c,"_vs_",s,".tsv")
Algo_output_file="DiamondX_vs_NR.tsv"
# DiamondXnr=fread(file =paste0(search_res_dir,"/","JGI_putative_RNA_Viruses_vs_ncbi_NR/DiamondX_JGI_putative_RNA_Viruses_vs_nr_eval0.01.tsv"), header = F, sep = "\t", col.names=Diamond_vNR_colnyms,stringsAsFactors = F)
DiamondXnr=fread(file=Algo_output_file, header = F, sep = "\t", col.names=Diamond_vNR_colnyms,stringsAsFactors = F)
# DiamondXnr=filter_by_LeP(Algo_out=DiamondXnr,min_len=min_len,pid=pid,evale=evale)
min_dt=data.frame(DiamondXnr[ , max(bits), by = qseqid])
colnames(min_dt)=c("qseqid","bits")
DiamondXnr=merge(min_dt,DiamondXnr,by=c("qseqid","bits"),all.x=T,all.y=F)
write.table(DiamondXnr,"Culled_DiamondX_vs_NR.tsv")
# DiamondXnr=DiamondXnr[-which(DiamondXnr$qseqid %in% FPs),]
dumas=grep(";",DiamondXnr$tax_id,fixed = T)

splitcoma <- function(x,y){
  return(unlist(strsplit(x["tax_id"],split=";", fixed=T))[y])
}

nodumas=DiamondXnr[-dumas,]
dumas_DiamondXnr=DiamondXnr[dumas,]
dumas_DiamondXnr$tax_id=apply(dumas_DiamondXnr, 1, splitcoma, 2)
DiamondXnr=rbind(nodumas,dumas_DiamondXnr)

merged_Diamonds=merge(DiamondXnr,lineages,by=c("tax_id"),all.x=T,all.y=F)
write.table(merged_Diamonds,"Culled_merged_Diamonds.tsv")

Cellular_Diamonds = merged_Diamonds[which(merged_Diamonds$superkingdom %in% Cellular_nms),]
provirus_relatives=c()
print("Started marking viral relatives")
for (i in Red_Herrings){ # Gluttony.
  viral_relatives_a=grep(i, Cellular_Diamonds[,c("salltitles")], ignore.case = T)
  viral_relatives_b=grep(i, Cellular_Diamonds[,c("sseqid")], ignore.case = T)
  provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
  print(paste0("Finished marking matches with the term: ",i))
}
# length(provirus_relatives)
if (isEmpty(provirus_relatives)==F){
  Cellular_Diamonds=Cellular_Diamonds[-provirus_relatives,]  
}
provirus_relatives=c()
for (i in Non_disqualifyiers){ # Gluttony.
  viral_relatives_a=grep(i, Cellular_Diamonds[,c("salltitles")], ignore.case = T)
  viral_relatives_b=grep(i, Cellular_Diamonds[,c("sseqid")], ignore.case = T)
  provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
  print(paste0("Finished marking matches with the term: ",i))
}
length(provirus_relatives)
if (isEmpty(provirus_relatives)==F){
  Cellular_Diamonds=Cellular_Diamonds[-provirus_relatives,]  
}
Cellular_Diamonds=filter_by_LeP(Algo_out=Cellular_Diamonds,min_len=(min_len/3),pid=(pid*100),evale=evale)
removies=unique(Cellular_Diamonds$qseqid)
keepers=setdiff(names(Input_fasta),removies)
output_fasta=Input_fasta[keepers]
writeXStringSet(output_fasta,"NR_filtered.fasta")

# FPs=unique(append(FPs,Cellular_Diamonds$qseqid))
# 
# ##### input vs NCBI's NT ##### 
# Algo="BLASTn"
# BLASTN_vNT=fread(file="cated_BLASTn_NT_vs_DiamondX_MG_filtered.tsv", header = F, sep = "\t", col.names=BLASTn_vNT_colnyms,stringsAsFactors = F)
# BLASTN_vNT=filter_by_LeP(Algo_out=BLASTN_vNT,min_len=min_len,pid=pid,evale=evale)
# BLASTN_vNT=BLASTN_vNT[-which(BLASTN_vNT$qseqid %in% FPs),]
# 
# dt <- data.table(BLASTN_vNT)
# min_dt=data.frame(dt[ , max(bits), by = qseqid]) #,
# colnames(min_dt)=c("qseqid","bits")
# BLASTN_vNT=merge(min_dt,BLASTN_vNT,by=c("qseqid","bits"),all.x=T,all.y=F)
# {
#   Cellular_Mtchs = BLASTN_vNT[which(BLASTN_vNT$sskingdoms %in% Cellular_nms),]
#   provirus_relatives=c()
#   print("Started marking viral relatives")
#   for (i in Red_Herrings){ # Gluttony.
#     viral_relatives_a=grep(i, Cellular_Mtchs[,c("stitle")], ignore.case = T)
#     viral_relatives_b=grep(i, Cellular_Mtchs[,c("sseqid")], ignore.case = T)
#     provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
#     print(paste0("Finished marking matches with the term: ",i))
#   }
#   if (isEmpty(provirus_relatives)==F){
#     Cellular_Mtchs=Cellular_Mtchs[-provirus_relatives,]
#   }
#   provirus_relatives=c()
#   for (i in Non_disqualifyiers){ # Gluttony.
#     viral_relatives_a=grep(i, Cellular_Mtchs[,c("stitle")], ignore.case = T,fixed = T)
#     viral_relatives_b=grep(i, Cellular_Mtchs[,c("sseqid")], ignore.case = T,fixed = T)
# 
#     provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
#     print(paste0("Finished marking matches with the term: ",i," which matched  ",(length(viral_relatives_b))))
#   }
#   if (isEmpty(provirus_relatives)==F){
#     Cellular_Mtchs=Cellular_Mtchs[-provirus_relatives,]
#   }
# }
# # FPs=unique(append(FPs,Cellular_Mtchs$qseqid))
# removies=unique(Cellular_Mtchs$qseqid)
# Input_fasta=readDNAStringSet("DiamondX_MG_filtered.fasta")
# keepers=setdiff(names(Input_fasta),removies)
# output_fasta=Input_fasta[keepers]
# writeXStringSet(output_fasta,"NR_filtered.fasta")

# ##### MMseqs2 vs Pfam32_A ##### 
evale=0.001
min_len=50
pid=0.25
MMseqs2_vPfam=fread("presence_absence/Search_results/MMseqs2_Pfam32_A_Results_Wdisc.tsv", header = T, sep = "\t", stringsAsFactors = F,data.table = F)
MMseqs2_vPfam=filter_by_LeP(Algo_out=MMseqs2_vPfam,min_len=min_len,pid=0.3,evale=evale)
# MMseqs2_vPfam=MMseqs2_vPfam[-which(MMseqs2_vPfam$qseqid %in% FPs),]
provirus_relatives=c()
print("Started marking viral relatives")
Red_Herrings=union(Red_Herrings,c("virid","RdRP"))
for (i in Red_Herrings){ # Gluttony.
  viral_relatives_a=grep(i, MMseqs2_vPfam[,c("PFAM_Name")], ignore.case = T)
  viral_relatives_b=grep(i, MMseqs2_vPfam[,c("PFAM_Long_verb")], ignore.case = T)
  viral_relatives_c=grep(i, MMseqs2_vPfam[,c("PFAM_desc")], ignore.case = T)
  
  provirus_relatives=append(provirus_relatives,union(viral_relatives_a,union(viral_relatives_c,viral_relatives_b)))
  print(paste0("Finished marking matches with the term: ",i))
}
if (isEmpty(provirus_relatives)==F){
  provirus_relatives_qseqids=unique(MMseqs2_vPfam$query[provirus_relatives])
  MMseqs2_vPfam=MMseqs2_vPfam[-which(as.character(MMseqs2_vPfam$query) %in% provirus_relatives_qseqids),]
}
repeat_relatives=c()
print("Started marking repeat relatives")
repeats_nms=c("repeat","Ank_","Ankyrins","PPR_","Fibronectin","PA14")
for (i in repeats_nms){ # Gluttony.
  repeat_relatives_a=grep(i, MMseqs2_vPfam[,c("PFAM_Name")], ignore.case = T)
  repeat_relatives_b=grep(i, MMseqs2_vPfam[,c("PFAM_Long_verb")], ignore.case = T)
  repeat_relatives_c=grep(i, MMseqs2_vPfam[,c("PFAM_desc")], ignore.case = T)
  
  repeat_relatives=append(repeat_relatives,union(repeat_relatives_a,union(repeat_relatives_b,repeat_relatives_c)))
  print(paste0("Finished marking matches with the term: ",i))
}
if (isEmpty(repeat_relatives)==F){
  repeat_df=MMseqs2_vPfam[repeat_relatives,]
}
hist_df=as.data.frame(MMseqs2_vPfam[,c("PFAM_Name","PFAM_Long_verb")] %>% group_by_all %>% count )
# repeat_df=hist_df[grep("repeat", hist_df[,c("PFAM_Long_verb")], ignore.case = T),]

Input_fasta=readDNAStringSet("presence_absence/Fasta_files/NT_filtered.fasta")

# dups=repeat_df$query[duplicated(repeat_df$query)]
# dous_repeat_df=repeat_df[which(repeat_df$query %in% dups),]
# removies=unique(dous_repeat_df$query)
removies=unique(repeat_df$query)

keepers=setdiff(names(Input_fasta),removies)
output_fasta=Input_fasta[keepers]
# writeXStringSet(output_fasta,"Repeat_filtered_1.fasta")
scaffold_df=data.frame("full_name"=names(output_fasta),"Length"=output_fasta@ranges@width,stringsAsFactors = F)
writeXStringSet(Input_fasta[union(removies,manual_FPs)],"removies_repeats.fasta")

##### Multirepeat-filtered vs repeatDB + New FPDB ##### 
min_len=90
pid=0.65
evale=0.000001
MMseqs2_vFPs=fread("./presence_absence/vsrepeatFPdb.tsv",col.names = MMseqs_vDNA_colnyms,stringsAsFactors = F)
MMseqs2_vFPs=filter_by_LeP(Algo_out=MMseqs2_vFPs,min_len=min_len,pid=pid,evale=evale)
min_dt=data.frame(MMseqs2_vFPs[ , max(bits), by = qseqid])
colnames(min_dt)=c("qseqid","bits")
MMseqs2_vFPs=merge(min_dt,MMseqs2_vFPs,by=c("qseqid","bits"),all.x=T,all.y=F)
# MMseqs2_vFPs=MMseqs2_vFPs[-which(duplicated(MMseqs2_vFPs$qseqid)),]
removies=unique(MMseqs2_vFPs$qseqid)
# removies=intersect(unique(MMseqs2_vFPs$qseqid),names(Input_fasta))
# Input_fasta=readDNAStringSet("Repeat_filtered_1.fasta")
Input_fasta=readDNAStringSet("/media/uri/HDD1/uri/V1_hard_mode/V2/V2.04/presence_absence/multirepeatFP2_fitered.fasta")
keepers=setdiff(names(Input_fasta),removies)

output_fasta=Input_fasta[keepers]
writeXStringSet(output_fasta,"multirepeatFP3_fitered.fasta")
scaffold_df=data.frame("full_name"=names(output_fasta),"Length"=output_fasta@ranges@width,stringsAsFactors = F)
MMseqs2_vPfam_cur=MMseqs2_vPfam[which(as.character(MMseqs2_vPfam$query) %in% scaffold_df$full_name),]
hist_df=as.data.frame(MMseqs2_vPfam_cur[,c("PFAM_Name","PFAM_Long_verb")] %>% group_by_all %>% count )
new_FPs=intersect(names(Input_fasta),union(manual_FPs,removies))
new_FPs2=setdiff(Input_fasta,output_fasta)

writeXStringSet(new_FPs2,"new_FPs3.fasta")

# ##### Re-examine the JGI scaf df ##### 
# filtered_DF=JGI_scaf_df[-which(JGI_scaf_df$full_name%in% FPs),]
# FPscafs=JGI_contigs[FPs]
# FPscafs1=intersect(JGI_contigs,FPscafs)
# 
# JGI_contigs_filt=setdiff(JGI_contigs,union(FPscafs1,reverseComplement(FPscafs1)))
# filtered_scafs=setdiff(names(JGI_contigs),names(JGI_contigs_filt))
# 
# 
