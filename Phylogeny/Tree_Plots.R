# Information  --------------------------------------------------------------------
## Script name: Tree_Plots.R
## Email: uri.neri@gmail.com
##

# Note! 
# This script was includes code that was not necessarily ran in a chronological order.
# Some of the objects were created by other Rscripts in the Phylogentic analysis (such as leaf2tax and cont2leaf).
# As such, it is not advised to attempt running this code 'as is', rather, it is intended to be a disclosure rather than a re-runnable "production" code.
# That said, we would be happy to try and help you should you want to source, modify, reuse, rerun or fork in any way the code.
# (Although, you don't have to tell us etc, this is open source so go ahead and use what you can as much as you want!)
# HTH
# Neri

# Body  --------------------------------------------------------------------
# #Plots  --------------------------------------------------------------------
library(ggplot2)
library(ggtree) 
library(ggrepel)
# library(plotly)
# library(ggplotify)
library(ggtreeExtra)
library(ggimage)
# library(treeio)
# library(tidytree)
library(ggstar)
library(ggnewscale)

###### Make ggtrees ###### 
# ggWorkTreeCirc = ggtree(WorkTree, layout = "circular", yscale = "none", ladderize = T)
# saveRDS(ggWorkTreeCirc, "ggWorkTreeCirc.phylo.tre.RDS")
ggWorkTreeCirc = readRDS("ggWorkTreeCirc.phylo.tre.RDS")

# ggWorkTreeApe = ggtree(WorkTree, layout = "ape", yscale = "none", ladderize = T)
# saveRDS(ggWorkTreeApe, "ggWorkTreeApe.phylo.tre.RDS")
ggWorkTreeApe = readRDS("ggWorkTreeApe.phylo.tre.RDS")

# ggWorkTreeRect = ggtree(WorkTree, layout = "rectangular", yscale = "none", ladderize = T)
# saveRDS(ggWorkTreeRect, "ggWorkTreeRect.phylo.tre.RDS")
ggWorkTreeRect = readRDS("ggWorkTreeRect.phylo.tre.RDS")


ggWorkTreeCircUltra <- readRDS("ggWorkTreeCircUltra.RDS.2022-01-12.17-25-13.RDS")
###### Prepare node annotation data ###### 

IDFT$BinID <- p0(
  str_split_fixed(string = IDFT$ND, pattern = fixed("_"),n = 3)[,1],
  "_",
  str_split_fixed(string = IDFT$ND, pattern = fixed("_"),n = 3)[,2])
IDFT[,MlenBin := sum(.SD),.SDcols="Length",by = (BinID)]


tmpAllDFs <- distinct(data.frame(rbindlist(list(AllDF,AllDF6Frx),fill=T))[,SharedCols(AllDF,AllDF6Frx)])
# Movement protein info:
w1 = grep(pattern = "MP_",x = tmpAllDFs$New_Name,fixed=T)
Putative_Movement_Encoding = unique(tmpAllDFs$ND[w1])
Putative_Movement_Encoding_Nodes = unique(na.omit(IDFT$node[which(IDFT$ND %in% Putative_Movement_Encoding)]))
w1 = intersect(wna(IDFT$Host),which(IDFT$node %in% Putative_Movement_Encoding_Nodes))
IDFT$Host[w1] <- "Eukaryota"
IDFT$Host_Evidence[w1] <- "Host related domain - Movement protein"


# Lysis info:
w1 = grep(pattern = "Lysis_sgl|Lysin_|LYZF1",x = tmpAllDFs$New_Name)
w2 = filter(tmpAllDFs[w1,], evalue < 0.005)
Putative_Lysis_Encoding = unique(w2$ND)
Putative_Lysis_Encoding_Nodes = unique(na.omit(IDFT$node[(intersect(wh(IDFT$ND %in% Putative_Lysis_Encoding),c(wh(IDFT$Phylum %in% c("Lenarviricota","p.0002.base-Kitrino")),wh(IDFT$Class == "Duplopiviricetes"))))]))

Putative_Lysis_Encoding_RID = unique(na.omit(IDFT$node.label[which(IDFT$ND %in% Putative_Lysis_Encoding)]))
w1 = intersect(wna(IDFT$Host),which(IDFT$node %in% Putative_Lysis_Encoding_Nodes))
IDFT$Host[w1] <- "Bacteria"
IDFT$Host_Evidence[w1] <- "Putative Host related domain - bacteriolytic protein"

# CRISPR info:
w1 = union(wana(IDFT$Hit.s.),wana(IDFT$Type.hit))
Putative_CRISPR_Matching_Nodes = unique(na.omit(IDFT$node[w1]))
IDFT$Host_Evidence[wh(IDFT$Host_Evidence == "Putative CRISPR spacer match")] = NA
IDFT$Host[wh(IDFT$Host_Evidence == "Putative CRISPR spacer match")] = NA
w1 = intersect(wna(IDFT$Host_Evidence),wh(IDFT$node %in% Putative_CRISPR_Matching_Nodes))
IDFT$Host[w1] <- "Bacteria"
IDFT$Host_Evidence[w1] <- "Putative CRISPR spacer match"

tmpIDFT <- IDFT
tmpIDFT$Novel <- T
tmpIDFT$Novel[grep(pattern = "al2_",x = tmpIDFT$RID,fixed = T)] = F
tmpIDFT$Novel[grep(pattern = "ap2_",x = tmpIDFT$RID,fixed = T)] = F
tmpIDFT$Novel[grep(pattern = "rv20_",x = tmpIDFT$RID,fixed = T)] <-   F
tmpIDFT$Novel[grep(pattern = "ya2_",x = tmpIDFT$RID,fixed = T)] = F


info <- distinct(tmpIDFT[wana(tmpIDFT$node),] %>% group_by(node) %>% summarise(
  node = node,
  Mlen = max(MlenBin),
  CRHit = toString(unique(na.omit(Type.hit))),
  Host = toString(unique(Host)),
  RBS = toString((RBS)),
  Phylum = unique(Phylum),
  Class = unique(Class),
  Order = unique(Order),
  Family = unique(Family),
  Genus = toString(unique(Genus)),
  # Novel = ((len(wh(Novel))) / len(Novel)  ),
  Novel = AllTrue(Novel),
  Genetic_Codes = toString(unique(Genetic_Code))))

info$node.label <- unlist(NID2NLabel(WorkTree,TmpNodeTblxIDDT$node))

info$Note <- NA
info$Note[info$node %in% Putative_Lysis_Encoding_Nodes]  <-  LysisSign
info$Note[info$node %in% Putative_CRISPR_Matching_Nodes] <-  CRISPRSign

saveRDS(info,"../Wolf/ggtreeinfo.12012022.RDS")

setDF(info)

infoInternals <- data.frame("Taxa"="","rank"="","node"="","node.label"="","maxlen"=0)[0,]
for (rnkcl in c("Phylum","Class","Order","Family","Genus")) {
  AllRnkTaxs <- data.frame(Taxa=unname(unique(info[,rnkcl])),rank=rnkcl)
  infoInternals <- rbind.fill(infoInternals,AllRnkTaxs)
}

setDT(infoInternals)

for (ix in 1:nrow(infoInternals)) {
  tips_ix <- unique(info$node[wh(info[,as.char(infoInternals$rank[ix])] == infoInternals$Taxa[ix])])
  tips_iy <- unlist(NID2NLabel(WorkTree,tips_ix))
  Mlen <- max(IDFT$MlenBin[wh(IDFT$node.label %in% tips_iy)])
  Taxlca <- get_mrca_of_set(WorkTree,tips_ix)
  infoInternals[ix,node:=Taxlca]
  infoInternals[ix,maxlen:=  Mlen]
}
infoInternals[,node := as.num(node)]
infoInternals[,node.label := unlist(NID2NLabel(WorkTree,node))]

infoInternals2 <- infoInternals[wh(infoInternals$rank != "Genus"),]
infoInternals3 <- infoInternals[grep(pattern = "genPartiti",x = infoInternals$Taxa),]
infoInternals <- rbind(infoInternals2,infoInternals3)
rm(infoInternals2,infoInternals3)

RTsLCA <- get_mrca_of_set(WorkTree,grep(x = WorkTree$tip.label,pattern = "rt",value = T))
infoInternals <- rbind(infoInternals,data.frame("node" = RTsLCA, "node.label" =NID2NLabel(WorkTree,RTsLCA),"rank"="Phylum","Taxa"="Reverse_Transcriptases","maxlen"=0))
w1 <- intersect(wh(infoInternals$rank != "Class"),grep(pattern = fixed(".base-"),x =infoInternals$Taxa,fixed = T))
infoInternals$NymN <- infoInternals$Taxa
infoInternals$NymN[w1] <- str_split_fixed(string = infoInternals$Taxa[w1],pattern = fixed(".base-"),n=2)[,1]
ToKeep <- intersect(wh(info$Note == "X"),unique(c(wh(info$Phylum %in% c("Lenarviricota","p.0002.base-Kitrino")),wh(info$Class == "Duplopiviricetes"))))

Nodetbl <- Tree2Edgetbl(WorkTree)
info2 <- rbind.fill(data.frame(info),data.frame(filter(Nodetbl,!(node.label %in% info$node.label)))[,c("node.label","node")])
infoInternals <- filter(infoInternals,!(Taxa  %in% c("Alvernaviridae", "Barnaviridae", "Carmotetraviridae", "Quadriviridae", "Sunviridae")))

setDT(info2)
for(ix in 1:nrow(infoInternals)){
  Allklids <- unlist(GetKids(WorkTree,node = as.num(infoInternals$node[ix]),type="all",Return_Labels = T))
  iy = wh(info2$node.label %in% Allklids)
  info2[iy, as.char(infoInternals$rank[ix]) :=  infoInternals$Taxa[ix]]
}
setDF(info2)
for(iy in 1:ncol(info2)){
  w1 <- wh(unlist(unname(info2[,iy])) == "NA")
  if(len(w1) !=0){
    info2[w1,iy] = NA
  }
}

infoInternals$DistFromRoot = get_all_distances_to_root(WorkTree, as_edge_count = F)[NLabel2NID(WorkTree, infoInternals$node.label)]
max(infoInternals$DistFromRoot)

info2$DistFromRoot = get_all_distances_to_root(WorkTree, as_edge_count = F)[NLabel2NID(WorkTree, info2$node.label)]
wna(info2$Novel)

AltSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/Alt_Sign.png")
PermuteSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/Permute_Sign.png")
CRISPRSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/BlueStar.png")
LysisSign <- p0("/media/HDD1/uri/RNA_Vir_MTs/misc/RedStar.png")

AltSets <- fread("Alt_Clades.tab",sep = "\t")
AltSets$leaflist = vector('list',c(nrow(AltSets)))   # List place holder 1
AltSets$leaflist <- apply(AltSets,MARGIN = 1,FUN = function(x) as.char(unlist(str_split_fixed(string = x["V2"],pattern = ",",n = Inf))))
AltSets$node <- unlist(apply(AltSets,MARGIN = 1,FUN = function(x) get_mrca_of_set(WorkTree,unlist(x["leaflist"]))))
AltSets$node.label <- NID2NLabel(WorkTree,AltSets$node)

PermutSets <- fread("Perm_Clades.txt",sep = "\t")
PermutSets$leaflist = vector('list',c(nrow(PermutSets)))   # List place holder 1
PermutSets$leaflist <- apply(PermutSets,MARGIN = 1,FUN = function(x) as.char(unlist(str_split_fixed(string = x["V2"],pattern = ",",n = Inf))))
PermutSets$node <- unlist(apply(PermutSets,MARGIN = 1,FUN = function(x) get_mrca_of_set(WorkTree,unlist(x["leaflist"]))))
PermutSets$node.label <- NID2NLabel(WorkTree,PermutSets$node)

info2$shapo <- NA
info2$shapo[wh(info2$node.label %in% AltSets$node.label)] <- "Alternative_Genetic_Code"
info2$shapo[wh(info2$node.label %in% PermutSets$node.label)] <- "Permuted_RdRP"

info2$Novel[wh(info2$Phylum == "Reverse_Transcriptases")] <- F
info2$Mlen[wh(info2$Phylum == "Reverse_Transcriptases")] <- 0


setDT(info2)
info2$AnyKnown <- !(info2$Novel)
info2$AnyKnown[wh(info2$Phylum == "Reverse_Transcriptases")] <- T
info2$NovelProp <- 0
info2$Total_weight <- 0
info2$Novel_weight = 0

info2$AllKids = GetKids(WorkTree,info2$node,type="tips")
info2$Nkids <- lengths(info2$AllKids)
info2$NkidsProp <- info2$Nkids/Ntip(WorkTree)

for(ix in 1:nrow(info2)){
  iy <- wh(info2$node %in% unlist(info2$AllKids[ix]))
  info2[ix, Total_weight:= sum(info2$DistFromRoot[iy])]
  info2[ix, AnyKnown:= AnyFalse(info2$Novel[iy])]
  iy <- wh(info2$Novel[iy])
  if(len(iy)==0){
    next
  }
  info2[ix,Novel_weight:= sum(info2$DistFromRoot[iy])]
  print(ix)
}

info2$NovelProp = info2$Novel_weight / info2$Total_weight
info2 <- distinct(info2)
saveRDS(info2,"../Wolf/ggtreeinfo.12012022.RDS")
# 
# WorkTree2$node.label = WorkTree$node.label
# ggWorkTreeCircUltra <- ggtree(WorkTree2,layout = "circular",  branch.length = "none")


###### MegaTree in Circular layout  - Color by novelty, label by taxonomic rank ###### 
library(ggtree)
library(ragg)

ggWorkTree = readRDS("ggWorkTreeCircUltra.RDS.2022-01-12.17-25-13.RDS")
ggWorkTree = readRDS("ggWorkTreeCirc.phylo.tre.RDS")

colorcols <- c("Lenarviricota" = "#0072B2", "Lenar" = "#0072B2",
               "Pisuviricota" = "#009E73","Pisu" = "#009E73",
               "Kitrinoviricota" = "#F0E442","Kitrino" = "#F0E442",
               "Duplornaviricota" = "#8E8C62","Duplorna" = "#8E8C62",
               "Negarnaviricota" = "#CC79A7","Negarna" = "#CC79A7",
               "p.0002.base-Kitrino" = "#DC143C","p.0002" = "#DC143C","<-----------p.0002----|" = "#DC143C",
               "Phylum" = "#DA8D7A","Family" = "#DAD97A", "Genus" = "#7ADA89",  "Order" = "#7ABADA",  "Class"= "#A17ADA",
               "Alt" = "Red", "Standard" = "Blue",
               "Alternative_Genetic_Code" = 17, "Permuted_RdRP"=16)


p2 <- ggWorkTree
p4 <- p2
p4$data$label <- unlist(NID2NLabel(WorkTree,p4$data$node))
p4$data <- merge((p4$data),info2,by = "node",all.x = T,all.y = F)
p4$data$AllKids <- NULL
p4$data <- merge((p4$data),`dropcols<-`(infoInternals,c("node.label","DistFromRoot","Lysis","CRISPR","Nkids","NkidsProp")),by = "node",all.x = T,all.y = F)
p4$data <- distinct(p4$data)
p4$data$Novel <- as.num(!(p4$data$AnyKnown))

## Color branches by novelty: 
p5 <- p4 + geom_tree(mapping = aes(color = as.num(AnyKnown))) + scale_color_gradient(low="#222222", high="#759495")# +

## Add Lysis/CRISPR tip labels (stars):
p6 <- p5  + scale_fill_manual(values = colorcols) +
  geom_fruit(geom=geom_image,mapping=aes(image = Note),
             pwidth=0.5,
             size = 0.0078,
             offset = 0.038) +
  theme(axis.title = element_text(hjust = -5,vjust = -4))

## Add Permuted RdRps and Alt. genetic codes arcs and text:
p7 <- p6  + geom_cladelab(
  mapping = aes(subset = (shapo=="Alternative_Genetic_Code"),
                node = node, label = RBS),
  barcolour = "#00FF11", offset = 0.34,label = "",barsize = 5,align = TRUE,fontsize=0, show.legend = T)

p8 <- p7  + geom_cladelab(
  mapping = aes(subset = (shapo=="Permuted_RdRP"),
                node = node, label = RBS),
  barcolour = "#EB8A23", offset = 0.56,label = "",barsize = 5,align = TRUE,fontsize=0, show.legend = T)

## Add Outer ring with barplot of genome lengths:
p9 <- p8 + scale_fill_manual(values = colorcols) +
  geom_fruit(geom=geom_bar,mapping=aes(x=Mlen,fill=Phylum),
             pwidth=0.25,orientation="y",stat="identity",
             grid.params=list(vline =F),
             axis.params=list(title = "Maximal observed Genome Length [nt]",
                              axis       = "x",
                              title.height = 0.001,
                              # text.size = 3.1,
                              hjust      = -0.2,
                              vjust      = +0.1,
                              nbreak     = 10,
                              line.size = 0.4),
             offset = 0.27) +
  theme(axis.title = element_text(hjust = -5,vjust = -4))
p9[["layers"]][[16]][["data"]][["y"]] <- 190

## Add Phyla arcs and text:
w1 <- Rename1Col(Rename1Col(filter(p9$data,rank == "Phylum"),"Taxa","cladelab"),"node","nodex")
w1$cladelab[wh(w1$cladelab=="p.0002.base-Kitrino")] <- "p.0002"
p10 <- p9 +#  scale_fill_manual(values = colorcols) + 
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Lenarviricota")],label =w1$cladelab[wh(w1$cladelab=="Lenarviricota")], barcolour="White",textcolour = unname(colorcols["Lenarviricota"]), offset = 0.728, align = T, offset.text = .02,hjust = "center", barsize = 0.0000000009, fontsize = 12, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Pisuviricota")], label =w1$cladelab[wh(w1$cladelab=="Pisuviricota")],barcolour="White",textcolour = unname(colorcols["Pisuviricota"]), offset = 0.708, align = T, offset.text = .02,hjust = "center", barsize = 0.000000010, fontsize = 12, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Kitrinoviricota")], label =w1$cladelab[wh(w1$cladelab=="Kitrinoviricota")],barcolour="White",textcolour = unname(colorcols["Kitrinoviricota"]), align = T,offset = 0.75, offset.text = .02,hjust = "center", barsize = 0.000000010, fontsize = 12, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Duplornaviricota")], label =w1$cladelab[wh(w1$cladelab=="Duplornaviricota")],barcolour="White",textcolour = unname(colorcols["Duplornaviricota"]), align = T, offset = 0.740, offset.text = .02,hjust = "center", barsize = 0.000000005, fontsize = 4.5, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Negarnaviricota")], label =w1$cladelab[wh(w1$cladelab=="Negarnaviricota")],barcolour="White",textcolour = unname(colorcols["Negarnaviricota"]), align = T, offset = 0.688, offset.text = .02,hjust = "center", barsize = 0.000000008, fontsize = 7, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="p.0002")], label =w1$cladelab[wh(w1$cladelab=="p.0002")],barcolour=unname(colorcols["p.0002"]),textcolour = unname(colorcols["p.0002"]), align = T, offset = 0.550, offset.text = .1,hjust = "center", barsize = 0.0000000005, fontsize = 11, angle = 329, horizontal = F, fontface="italic") 


Groups2Tag <- p10$data$Taxa[wh(p10$data$rank == "Order")]
Groups2Tag <- Groups2Tag[-grep(pattern = fixed("."),x = Groups2Tag,fixed = T)]
Groups2Tag <- unique(c(Groups2Tag,"Reverse_Transcriptases","Hepelivirales","Pisuviricota","Duplornaviricota","p.0002.base-Kitrino","Negarnaviricota","f.0278","Lenarviricota","Kitrinoviricota","genPartiti.0019.base-Deltapartitivirus", "Cystoviridae","Picobirnaviridae","Partitiviridae","f.0271.base-Toga","f.0273.base-Toga","f.0268.base-Toga","f.0066.base-Hypo","f.0069.base-Hypo","f.0068.base-Hypo","f.0215.base-Deltaflexi","Nidovirales","Tobaniviridae","Hypoviridae","Deltaflexiviridae","Benyviridae","Closteroviridae","Botourmiaviridae","f.0024.base-Solinvi","Polycipiviridae","Rhabdoviridae","Coronaviridae"))
Groups2Tag <- setdiff(Groups2Tag,c("Muvirales"))
p10$data$NymN[wh(!(p10$data$Taxa %in% Groups2Tag))] = NA

p10$data$NymN <- gsub(x = p10$data$NymN,replacement = "",pattern = gsub(x = toString(RankDF$RegVec),pattern = ", ", replacement = "|"))
p10$data <- distinct(p10$data)
p11 <- p10 +  geom_label_repel(aes(label=NymN,fill=rank), segment.linetype = 6, segment.colour="Red", label.size = 0, max.time = 1, max.iter = 100, max.overlaps = 30, label.padding = 0.2, size = 5, na.rm = T,
                               fontface="italic",box.padding	=0,
                               arrow = arrow(length = unit(0.01, "cm"),type="closed",ends = "last"))

### Putting labels on the Phyla
w1 <- filter(p11$data,rank == "Phylum")
w1 <- Rename1Col(w1,"Taxa","cladelab")
w1 <- Rename1Col(w1,"node","nodex")
w1$cladelab[wh(w1$cladelab=="p.0002.base-Kitrino")] <- "<----p.0002----|"
p12 <- p11 +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Lenarviricota")],label =w1$cladelab[wh(w1$cladelab=="Lenarviricota")], barcolour=unname(colorcols["Lenarviricota"]), offset = -0.8, align = T, offset.text = .02,hjust = "center", barsize = 0.9, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Pisuviricota")], label =w1$cladelab[wh(w1$cladelab=="Pisuviricota")],barcolour=unname(colorcols["Pisuviricota"]), offset = -1.1, align = T, offset.text = .02,hjust = "center", barsize = 1.0, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Kitrinoviricota")], label =w1$cladelab[wh(w1$cladelab=="Kitrinoviricota")],barcolour=unname(colorcols["Kitrinoviricota"]), align = T,offset = -0.5, offset.text = .02,hjust = "center", barsize = 1.0, fontsize = 9, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Duplornaviricota")], label =w1$cladelab[wh(w1$cladelab=="Duplornaviricota")],barcolour=unname(colorcols["Duplornaviricota"]), align = T, offset = -1, offset.text = .02,hjust = "center", barsize = .5, fontsize = 4, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="Negarnaviricota")], label =w1$cladelab[wh(w1$cladelab=="Negarnaviricota")],barcolour=unname(colorcols["Negarnaviricota"]), align = T, offset = -1.2, offset.text = .02,hjust = "center", barsize = .8, fontsize = 6, angle = "auto", horizontal = F, fontface="italic") +
  geom_cladelab(node = w1$nodex[wh(w1$cladelab=="<----p.0002----|")], label =w1$cladelab[wh(w1$cladelab=="<----p.0002----|")],barcolour=unname(colorcols["p.0002"]), align = T, offset = -1, offset.text = .1,hjust = "center", barsize = .5, fontsize = 6, angle = 330, horizontal = F, fontface="italic") 


p12$data$Alt=NA
p12$data$Alt[wh(p12$data$node %in% AltSets$node)] <- T

p12$data$Perm=NA
p12$data$Perm[wh(p12$data$node %in% PermutSets$node)] <- T

p12$data$shapo = NA
p12$data$shapo[wh(p12$data$Alt)] <- "Alternative Genetic Code"
p12$data$shapo[wh(p12$data$Perm)] <- "Permuted RdRP"

p12$data$NkidsProp[] <-(16+3*(scales::rescale(x = as.integer(p12$data$Nkids),to=c(1,2))))

p13 <- p12 + geom_point2(aes(size=(6+(NkidsProp*24)),shape=shapo), na.rm = T,color = "purple", fill='purple')
p13 <- p12 +   geom_nodepoint(mapping = aes(subset = (!is.na(shapo)),shape=shapo,size = (8+(NkidsProp*24)),fill=shapo),)#,color = "purple", fill='purple') #+




# tmppppp <- p13
# tmppppp$layers[[1]] = NULL
# tmppppp$data <- `dropcols<-`(tmppppp$data,c("DistFromRoot.x","Lysis.y","Lysis.x","CRISPR.x","CRISPR.y","DistFromRoot.y"))
# fastSave::saveRDS.lbzip2(n.cores = 12,object = tmppppp,file = p0("AlignedTip.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".RDS"))
# {
#   source("/media/HDD1/uri/RNA_Vir_MTs/RVMT/Misc/basicf.r")
#   setwd("/media/HDD1/uri/RNA_Vir_MTs/V3/Wolf/workdir/")
#   library(ggtreeExtra);library(ggstar);library(ggnewscale);library(ggtree);library(stringr);library(ggimage);library(gginnards)
#   rm(testasd)
#   testasd <- QuickRDSRead(n.cores=11,fname="tmp.*.RDS")
#   testasd <- fastSave::readRDS.lbzip2(n.cores = 12,file = list.files(path = ".",pattern = "^AlignedTip.*RDS"))
#   pdf(encoding = "ISOLatin1.enc",file = paste0("AlignedTip.228.",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".QuickTree.pdf"),width = 41, height = 41,title = "Riboviria (CRISPR+Lysis stared, Permutation + Alt code as internal shapes)",useDingbats = T)
#   testasd# +theme(legend.position="bottom")
#   dev.off()
#   
# }
