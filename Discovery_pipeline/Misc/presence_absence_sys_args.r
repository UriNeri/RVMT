# Presence Absence (generic).
##### Set env ##### 
library(Biostrings)
library(readr)
library(data.table)
source(functions_file)

##### Arguments  ##### 
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
  Non_disqualifyiers_file="/urigo/urineri/RNA_Viruses/acc/ND_list.txt"
  functions_file="/urigo/urineri/scripts/RNA_Viruses/v2/R_scripts/fun_ctions.r"
  lineages <- fread("/urigo/urineri/RNA_Viruses/acc/lineages_merged.tsv") # Get all for now, later maybe keep only Kingdome
} # Args.
{ # Args for interactive jobs.
  Algo="DiamondX"  # BlastN, DiamondX, or DC-MegaBlast.
  iter=7 # Iteration
  pid=45 #Precent iden. 
  evale="0.1" #E-value
  scov=-1 # Subject coverage. Set to a negative value if you didn't set it for the Algo run (which opts to DEF).
  qcov=-1  # Query coverage. Set to a negative value if you didn't set it for the Algo run (which opts to DEF).
  # PWD="/scratch200/urineri/V2/absenties/"   # Current working directory/outdir. Something like "/scratch200/urineri/V2/output/" Keep the leading and trailing "/" !!!
  min_len=30# Minimal alignment length.
  DB="NR" # Subject database.
  min_bitscore=82 # only for DC-MegaBlast.
  Non_disqualifyiers_file="/media/uri/HDD1/uri/V1_hard_mode/V2/acc/ND_list.txt"
  functions_file="/urigo/urineri/scripts/RNA_Viruses/v2/R_scripts/fun_ctions.r"
  PWD="/media/uri/HDD1/uri/V1_hard_mode/"
  lineages <- fread(paste0(PWD,"/V2/acc/lineages_merged.tsv"),na.strings = "NA") # Get all for now, later maybe keep only Kingdome
}
##### Set consts.  ##### 
Non_disqualifyiers=scan(Non_disqualifyiers_file, what="", sep="\n",)
Cellular_nms= c("Archaea","Eukaryota","Bacteria")
Red_Herrings=c("virus","phage","capsid","coat","tail","LysM","lysis", "viral","vector","AngRem","PRSV","ds-RNA element","dsRNA element")
constructs_taxIDs=c("32630")
Diamond_vMGs_colnyms=c("qseqid", "sseqid", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "length", "pident", "mismatch", "slen")
BLASTn_vMGs_colnyms=c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore")
BLASTn_vNT_colnyms=c("qseqid", "sseqid", "sscinames", "pident", "qlen", "length", "mismatch", "evalue", "bitscore", "stitle", "sskingdoms", "tax_id")
Diamond_vNR_colnyms=c("qseqid", "sseqid", "tax_id", "stitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue","bitscore","salltitles")

##### Get to it  ##### 

setwd(PWD)
##### Presence - Absence ##### 
if (DB=="MGs"){
  study_list=(list.dirs(full.names = F,recursive = F))
  for (study in study_list) {
    print(paste0("Started  ",study))
    if(Algo=="DiamondX"){
      colnyms=Diamond_vMGs_colnyms
      if (qcov<0){
        output_file=paste0(PWD,study,"/",Algo,"_qMT_vs_sMG_scov",scov,"_eval",evale,".tsv")
      }
      if (scov<0){
        output_file=paste0(PWD,study,"/",Algo,"_qMT_vs_sMG","_qcov",qcov,"_eval",evale,".tsv")
      }
    }
    if(Algo=="BlastN"){
      colnyms=BLASTn_vMGs_colnyms
      output_file=paste0(PWD,study,"/",Algo,"_qMT_vs_sMG_Pid",pid,"_qcov",qcov,"_eval",evale,".tsv")
    }
    if (iter==0){
      scaffolds_fasta = readDNAStringSet(file = paste0("/scratch200/urineri/V2/input/MTs/",study,"/rmdup_cated_qualified_scaffolds_L>1200GN=>0.fasta"))
    }
    if (iter!=0){
        scaffolds_fasta = readDNAStringSet(file = paste0(PWD,study,"/absent_scafs_",study,"_iter",iter,".fasta"))
    }
    Algo_out <- try(fread(file =output_file, header = F, sep = "\t", col.names=colnyms))
    if (class(Algo_out)!="try-error"){
        meta_scaf=data.frame("scaf_id"=names(scaffolds_fasta),stringsAsFactors = F)
        Discarded_scafs=as.character(unique(Algo_out$qseqid[which(Algo_out$length>min_len)]))
        all_scafs=as.character(unique(meta_scaf$scaf_id))
        absent_scafs=scaffolds_fasta[setdiff(all_scafs,Discarded_scafs)]
        writeXStringSet(absent_scafs, paste0(PWD,study,"/absent_scafs_",study,"_iter",(iter+1),".fasta"))
        print(paste0("For ", study,"  the pass/total ratio is:   ",(length(absent_scafs)/(length(meta_scaf$scaf_id))),", (passed =  ",length(absent_scafs),")"))
      }
      if (class(Algo_out)=="try-error"){
        print("Algo_out file is empty, absent_scafs_blastndiamond_ would be the same as absent_study")
        writeXStringSet(scaffolds_fasta, paste0(PWD,study,"/absent_scafs_",study,"_iter",(iter+1),".fasta"))
      }
  }
} 

##### Refrence based clean - up ##### 
if (DB!="MGs"){
  if (iter==0){
    scaffolds_fasta = readDNAStringSet(file = paste0("/scratch200/urineri/V2/input/MTs/",study,"/rmdup_cated_qualified_scaffolds_L>1200GN=>0.fasta"))
  }
  if (iter!=0){
    scaffolds_fasta = readDNAStringSet(file = paste0(PWD,"cated_absenties_iter",iter,".fasta"))
  }
  print(paste0("Started fitlering absenties by ",Algo," vs ",DB))
  if(Algo=="DiamondX"){
    if (qcov<0){
      output_file=paste0(PWD,Algo,"_cated_absenties_iter",iter,"_vs_",DB,"_scov",scov,"_eval",evale,".tsv")
    }
    if (scov<0){
      output_file=paste0(PWD,Algo,"_cated_absenties_iter",iter,"_vs_",DB,"_qcov",qcov,"_eval",evale,".tsv")
    }
    colnyms=Diamond_vNR_colnyms
    Algo_out <- try(fread(file=output_file, header = F, sep = "\t", col.names=colnyms,stringsAsFactors = F))
    if (class(Algo_out)!="try-error"){
      dumas=grep(";",Algo_out$tax_id,fixed = T)
      for (i in dumas){ 
        Algo_out$tax_id[i]=unlist(strsplit(Algo_out$tax_id[i],split=";", fixed=T))[2]
      }
      Algo_out$tax_id=as.numeric(Algo_out$tax_id)
      print(paste0("Finished numerifying the input file"))
      merged_Diamonds=merge(Algo_out,lineages,by=c("tax_id"),all.x=T,all.y=F)
      print(paste0("Finished merging the Algo_out file with the lineages file."))
      {
        Cellular_Diamonds = merged_Diamonds[which(merged_Diamonds$superkingdom %in% Cellular_nms),]
        provirus_relatives=c()
        print("Started marking viral relatives")
        for (i in Red_Herrings){ # Gluttony.
          viral_relatives_a=grep(i, Cellular_Diamonds[,c("salltitles")], ignore.case = T)
          viral_relatives_b=grep(i, Cellular_Diamonds[,c("sseqid")], ignore.case = T)
          provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
          print(paste0("Finished marking matches with the term: ",i))
        }
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
        if (isEmpty(provirus_relatives)==F){
          Cellular_Diamonds=Cellular_Diamonds[-provirus_relatives,]  
        }
        Cellular_Diamonds=Cellular_Diamonds[which(Cellular_Diamonds$length>min_len),]
        Cellular_Diamonds=Cellular_Diamonds[which(Cellular_Diamonds$pident>pid),]
        Discarded_scafs=as.character(Cellular_Diamonds$qseqid)
      }
    }
    if (class(Algo_out)=="try-error"){
      print("Algo_out file is empty, resulting fasta would be the same as input.")
      writeXStringSet(scaffolds_fasta, paste0(PWD,study,"/absent_scafs_",study,"_iter",(iter+1),".fasta"))
    }
   }
  if(Algo %in% c("BlastN","DC-MegaBlast")){
    colnyms=BLASTn_vNT_colnyms
    if(Algo =="BlastN"){
      output_file=paste0(PWD,Algo,"_cated_absenties_iter",iter,"_vs_",DB,"_Pid",pid,"_scov",scov,"_eval",evale,".tsv")
    }
    if(Algo =="DC-MegaBlast"){
      output_file=paste0(PWD,Algo,"_cated_absenties_iter",iter,"_vs_",DB,"_eval",evale,".tsv")
    }
    Algo_out <- try(fread(file =output_file, header = F, sep = "\t", col.names=colnyms,stringsAsFactors = F))
    if (class(Algo_out)!="try-error"){
      if(Algo =="DC-MegaBlast"){ #Culling
        dt <- data.table(Algo_out)
        min_dt=data.frame(dt[ , max(bitscore), by = qseqid]) #,
        colnames(min_dt)=c("qseqid","bitscore")
        Algo_out=merge(min_dt,Algo_out,by=c("qseqid","bitscore"),all.x=T,all.y=F)
      }
      {
        Cellular_Mtchs = Algo_out[which(Algo_out$sskingdoms %in% Cellular_nms),]
        provirus_relatives=c()
        print("Started marking viral relatives")
        for (i in Red_Herrings){ # Gluttony.
          viral_relatives_a=grep(i, Cellular_Mtchs[,c("stitle")], ignore.case = T)
          viral_relatives_b=grep(i, Cellular_Mtchs[,c("sseqid")], ignore.case = T)
          provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
          print(paste0("Finished marking matches with the term: ",i))
        }
        if (isEmpty(provirus_relatives)==F){
          Cellular_Mtchs=Cellular_Mtchs[-provirus_relatives,]  
        }
        provirus_relatives=c()
        for (i in Non_disqualifyiers){ # Gluttony.
          viral_relatives_a=grep(i, Cellular_Mtchs[,c("stitle")], ignore.case = T)
          viral_relatives_b=grep(i, Cellular_Mtchs[,c("sseqid")], ignore.case = T)
          
          provirus_relatives=append(provirus_relatives,union(viral_relatives_a,viral_relatives_b))
          print(paste0("Finished marking matches with the term: ",i))
        }
        if (isEmpty(provirus_relatives)==F){
          Cellular_Mtchs=Cellular_Mtchs[-provirus_relatives,]  
        }
        Cellular_Mtchs=Cellular_Mtchs[which(Cellular_Mtchs$length>min_len),]
        Cellular_Mtchs=Cellular_Mtchs[which(Cellular_Mtchs$pident>pid),]
        Cellular_Mtchs=Cellular_Mtchs[which(Cellular_Mtchs$bitscore>min_bitscore),]
        Discarded_scafs=as.character(Cellular_Mtchs$qseqid)
      }
    }
    if (class(Algo_out)=="try-error"){
      print("Algo_out file is empty, resulting fasta would be the same as input.")
      Discarded_scafs=c("none_real")
    }
  }
  all_scafs=as.character(unique(names(scaffolds_fasta)))
  absent_scafs=scaffolds_fasta[setdiff(all_scafs,Discarded_scafs)]
  writeXStringSet(absent_scafs, paste0(PWD,"cated_absenties_iter",(iter+1),".fasta"))
  print(paste0("For iter:   ", iter,"  the pass/total ratio is:   ",(length(absent_scafs)/(length(scaffolds_fasta))),", (passed =  ",length(absent_scafs),")"))
  }

timestamp()


##### misc. ##### 

{  # Get matches from previous iterations.
#   scaf_data=data.frame("qseqid"=names(absent_scafs),"Length"=absent_scafs@ranges@width,stringsAsFactors = F)
#   scaf_data_1=merge(scaf_data,Algo_out,by=1,all.x=T,all.y=F)
#   Dcolnyms=c("qseqid", "sseqid", "staxids", "stitle", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue","bitscore","salltitles")
#   Bcolnyms=c("qseqid", "sseqid", "sscinames", "pident", "qlen", "length", "mismatch", "evalue", "bitscore", "stitle", "sskingdoms", "staxids")
#   Diamond_nr_cated=paste0(PWD,"cated_diamonds.tsv")
#   BlastN_nt_cated=paste0(PWD,"cated_blastn.tsv")
#   Diamond_nr_out <- try(fread(file=Diamond_nr_cated, header = F, sep = "\t", col.names=Dcolnyms,stringsAsFactors = F))
#   for (i in 1:nrow(Diamond_nr_out)){ # Heavy, but need make sure multiple hits won't get into NAs.
#     Diamond_nr_out$staxids[i]=as.numeric(unlist(strsplit(Diamond_nr_out$staxids[i],split=";", fixed=T))[1])
#   }
#   BlastN_nt_out <- try(fread(file=BlastN_nt_cated, header = F, sep = "\t", col.names=Bcolnyms,stringsAsFactors = F))
#   lineages_kingdonly=lineages[,c(1,2)]
#   colnames(lineages_kingdonly)=c("staxids","superkingdom")
#   Diamond_nr_out=merge(Diamond_nr_out,lineages_kingdonly,by=c("staxids"),all.x=T,all.y=F)
#   # rm(merged_Diamonds)
#   scaf_data=data.frame("old_ID"=names(absent_scafs),"Length"=absent_scafs@ranges@width,stringsAsFactors = F)
#   scaf_data_nr=merge(scaf_data,Diamond_nr_out,by=1,all.x=T,all.y=F)
#   scaf_data_nr_no_viruses=scaf_data_nr[-which(scaf_data_nr$ID %in% scaf_data_nr$ID[which(scaf_data_nr$s)] )]
#   scaf_data_nt=merge(scaf_data,BlastN_nt_out,by=1,all.x=T,all.y=F)
#   # Largest_scaf = absent_scafs[which(absent_scafs@ranges@width==max(absent_scafs@ranges@width))]
#   # writeXStringSet(Largest_scaf, paste0(PWD,"Largest_scaf",(iter+1),".fasta"))
#   # writeXStringSet(absent_scafs["3300012699_Ga0157593_1042799"] , paste0(PWD,"3300012699_Ga0157593_1042799.fasta"))
#  scaffolds_IDs_MTs_dffed <- read_delim("rmdup1200_scaffolds_IDs_MTs_dffed.txt","\t", escape_double = F, col_names = F)
#  scaffolds_IDs_MTs_dffed=scaffolds_IDs_MTs_dffed[,-3]
#  colnames(scaffolds_IDs_MTs_dffed)=c("new_ID","old_ID")
#  scaf_data=merge(scaf_data,scaffolds_IDs_MTs_dffed,by=c("old_ID"),all.x=T,all.y=F)
#  names(absent_scafs)=scaf_data$new_ID
#  writeXStringSet(absent_scafs, paste0(PWD,"V2_Scaffolds.fasta"))
}
                                       