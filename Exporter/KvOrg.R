# Information  --------------------------------------------------------------------
## Script name: KvOrg.R
## Email: uri.neri@gmail.com
# Body  --------------------------------------------------------------------
source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
RVMT <- "/media/HDD1/uri/RNA_Vir_MTs/RVMT/"
RNA_Vir_MTs_dir <- "/media/HDD1/uri/RNA_Vir_MTs/V3/"
versiza <- "Ben"
current_path <- p0(RNA_Vir_MTs_dir, versiza, "/")
setwd(current_path)
library(ape)
library(castor)
library(phangorn)




THREADS <- 11
Memory <- 11800


WriteScaffoldCart <- function(dfdt,output2){
  tmpdfdt <- dfdt[,c("Source","Full_name")]
  tmpdfdt$ScaffoldID <- str_split_fixed(string = tmpdfdt$Full_name,pattern = "_",n=2)[,-1]
  tmpdfdt$Full_name <- "assembled"
  colnames(tmpdfdt) <- c("Scaffold ID","","")
  cat("Scaffold ID",sep = "\n",file = output2)
  fwrite(x = tmpdfdt,file = output2,sep = " ",quote = F,append = T)
}


IDFT4Expo <- IDFT
IDFT4Expo$Phylum[wna(IDFT4Expo$Phylum)] <- "Unclassified"
IDFT4Expo <- `dropcols<-`(IDFT4Expo,c("Matt","RN90","NC90","NC95","NC98","NC99","DupNames","note","Was_In_Old_IDFT"))


tmpAllDFs4Expo <- filter(tmpAllDFs,ND %in% IDFT2$ND)
tmpAllDFs4Expo <- `dropcols<-`(tmpAllDFs,c("Family","Phylum","Genus","Class","Order","Classified","PolyPort_Problematic","Analysis_Type"))

WorkTree4Expo <- WorkTree
WorkTree4Expo$node.label <- NULL

setwd(current_path)
# unlink("RiboV1.4",recursive = T)
dir.create("RiboV1.4")
setwd("RiboV1.4")


# Let's go code duplication !

WriteWolfTbl(IDFT4Expo,p0("RiboV1.4_Info.tsv"))
WriteScaffoldCart(filter(IDFT4Expo,Novel == T), p0("RiboV1.4","_IMG_Scaffold_cart.tsv"))
writeXStringSet(IDFT.fasta[intersect(names(IDFT.fasta),IDFT4Expo$ND)],p0("RiboV1.4","_Contigs.fasta.gz"),compress = T)
WriteWolfTbl(filter(tmpAllDFs,ND%in% IDFT4Expo$ND),p0("RiboV1.4_HMMatches.tsv"))
try(write.tree(get_subtree_with_tips(WorkTree4Expo,only_tips = unique(IDFT4Expo$node.label))$subtree,file = p0("RiboV1.4","_subtree.newick")))
zip::zip(zipfile = p0("RiboV1.4",".zip"),files =c(p0("RiboV1.4","_subtree.newick"),p0("RiboV1.4","_HMMatches.tsv"),p0("RiboV1.4","_Contigs.fasta.gz"),p0("RiboV1.4","_Info.tsv"),p0("RiboV1.4","_IMG_Scaffold_cart.tsv")) )

for(Phyla in unique(na.omit(IDFT4Expo$Phylum))){
  # for(Phyla in setdiff(unique(na.omit(IDFT4Expo$Phylum)),c("Lenarviricota"))){

  print(Phyla)
  dir.create(Phyla)
  setwd(Phyla)
  PhyluDF <- filter(IDFT4Expo,Phylum == Phyla)
  
  WriteWolfTbl(PhyluDF,p0(Phyla,"_Info.tsv"))
  WriteScaffoldCart(filter(PhyluDF,Novel == T), p0(Phyla,"_IMG_Scaffold_cart.tsv"))
  writeXStringSet(IDFT.fasta[intersect(names(IDFT.fasta),PhyluDF$ND)],p0(Phyla,"_Contigs.fasta.gz"),compress = T)
  WriteWolfTbl(filter(tmpAllDFs,ND%in% PhyluDF$ND),p0(Phyla,"_HMMatches.tsv"))
  if(Phyla == "Unclassified"){
    zip::zip(zipfile = p0(Phyla,".zip"),files =c(p0(Phyla,"_HMMatches.tsv"),p0(Phyla,"_Contigs.fasta.gz"),p0(Phyla,"_Info.tsv"),p0(Phyla,"_IMG_Scaffold_cart.tsv")) )
    setwd("..")
    print(getwd())
    next
  }
  try(write.tree(get_subtree_with_tips(WorkTree4Expo,only_tips = unique(PhyluDF$RCR90))$subtree,file = p0(Phyla,"_subtree.newick")))
  zip::zip(zipfile = p0(Phyla,".zip"),files =c(p0(Phyla,"_subtree.newick"),p0(Phyla,"_HMMatches.tsv"),p0(Phyla,"_Contigs.fasta.gz"),p0(Phyla,"_Info.tsv"),p0(Phyla,"_IMG_Scaffold_cart.tsv")) )
  
  
  for(Clasa in unique(na.omit(PhyluDF$Class))){
    print(Clasa)
    dir.create(Clasa)
    setwd(Clasa)
    ClasaDF <- filter(PhyluDF,Class == Clasa)
    
    WriteWolfTbl(ClasaDF,p0(Clasa,"_Info.tsv"))
    WriteScaffoldCart(filter(ClasaDF,Novel == T), p0(Clasa,"_IMG_Scaffold_cart.tsv"))
    writeXStringSet(IDFT.fasta[ClasaDF$ND],p0(Clasa,"_Contigs.fasta.gz"),compress = T)
    WriteWolfTbl(filter(tmpAllDFs,ND%in% ClasaDF$ND),p0(Clasa,"_HMMatches.tsv"))
    x <-  try(write.tree(get_subtree_with_tips(WorkTree4Expo,only_tips = unique(ClasaDF$RCR90))$subtree,file = p0(Clasa,"_subtree.newick")))
    if (class(x) == "try-error"){
      system(command = p0('echo "Can not print subtree for singltons " >',p0(Clasa,"_subtree.newick")))
    }
    zip::zip(zipfile = p0(Clasa,".zip"),files =c(p0(Clasa,"_subtree.newick"),p0(Clasa,"_HMMatches.tsv"),p0(Clasa,"_Contigs.fasta.gz"),p0(Clasa,"_Info.tsv"),p0(Clasa,"_IMG_Scaffold_cart.tsv")) )
    
    for(Ordy in unique(na.omit(ClasaDF$Order))){
      print(Ordy)
      print(getwd())
      dir.create(Ordy)
      setwd(Ordy)
      OrdyDF <- filter(ClasaDF,Order == Ordy)
      
      WriteWolfTbl(OrdyDF,p0(Ordy,"_Info.tsv"))
      WriteScaffoldCart(filter(OrdyDF,Novel == T), p0(Ordy,"_IMG_Scaffold_cart.tsv"))
      writeXStringSet(IDFT.fasta[OrdyDF$ND],p0(Ordy,"_Contigs.fasta.gz"),compress = T)
      WriteWolfTbl(filter(tmpAllDFs,ND%in% OrdyDF$ND),p0(Ordy,"_HMMatches.tsv"))
      x <- try(write.tree(get_subtree_with_tips(WorkTree4Expo,only_tips = unique(OrdyDF$RCR90))$subtree,file = p0(Ordy,"_subtree.newick")))
      if (class(x) == "try-error"){
        system(command = p0('echo "Can not print subtree for singltons " >',p0(Ordy,"_subtree.newick")))
      }
      zip::zip(zipfile = p0(Ordy,".zip"),files =c(p0(Ordy,"_subtree.newick"),p0(Ordy,"_HMMatches.tsv"),p0(Ordy,"_Contigs.fasta.gz"),p0(Ordy,"_Info.tsv"),p0(Ordy,"_IMG_Scaffold_cart.tsv")) )
      
      for(Familia in unique(na.omit(OrdyDF$Family))){
        print(Familia)
        dir.create(Familia)
        setwd(Familia)
        FamiliaDF <- filter(OrdyDF,Family == Familia)
        WriteWolfTbl(FamiliaDF,p0(Familia,"_Info.tsv"))
        WriteScaffoldCart(filter(FamiliaDF,Novel == T), p0(Familia,"_IMG_Scaffold_cart.tsv"))
        writeXStringSet(IDFT.fasta[FamiliaDF$ND],p0(Familia,"_Contigs.fasta.gz"),compress = T)
        WriteWolfTbl(filter(tmpAllDFs,ND%in% FamiliaDF$ND),p0(Familia,"_HMMatches.tsv"))
        x <- try(write.tree(get_subtree_with_tips(WorkTree4Expo,only_tips = unique(FamiliaDF$RCR90))$subtree,file = p0(Familia,"_subtree.newick")))
        if (class(x) == "try-error"){
          system(command = p0('echo "Can not print subtree for singltons " >',p0(Familia,"_subtree.newick")))
        }
        try(zip::zip(zipfile = p0(Familia,".zip"),files =c(p0(Familia,"_subtree.newick"),p0(Familia,"_HMMatches.tsv"),p0(Familia,"_Contigs.fasta.gz"),p0(Familia,"_Info.tsv"),p0(Familia,"_IMG_Scaffold_cart.tsv")) ))
        setwd("..")
      }
      setwd("..")
    }
    setwd("..")
  }
  setwd("..")
  print(getwd())
}
  



