# Information  --------------------------------------------------------------------
## Script name: basicf.r
## Description: basic functions, stable across shenanigans... Contains a lot of bad code for legacy and nostaligic reasons only.
## Copyright (c) Uri Neri, 2020-2022
## Email: uri.neri@gmail.com

# set glob options  --------------------------------------------------------------------
options(scipen = 6, digits = 22) # I prefer to view outputs in non-scientific notation
options(stringsAsFactors = FALSE)
# on.exit(options(opts), add = TRUE);

# Install packages  --------------------------------------------------------------------
# WARNING! this next part is kind of invasive... tries to install packages from who knows where :-)
# # From CRAN:
# LibListx <- c("Rcpp", "castor", "phytools", "devtools", "methods", "plyr", "dplyr", "data.table", "stringr", "ggrepel", "ggplot2", "phangorn")
# for (lib in LibListx) {
#   x <- try(library(package = as.character(lib), character.only = T))
#   if (class(x) == "try-error") {
#     install.packages(pkgs = lib)
#   }
# }
# # From GitHub:
# devtools::install_github("jimhester/knitrBootstrap")
# devtools::install_github("barkasn/fastSave")
# devtools::install_github("urineri/biofiles")
# # devtools::install_github("thackl/thacklr")
# # devtools::install_github("thackl/gggenomes")
# devtools::install_github("berndbischl/BBmisc")


# # From BioConducter:
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages(pkgs = "BiocManager")
# }

# LibListx <- c("Biostrings", "DECIPHER", "ggtree", "GenomicRanges", "ggtree", "plyranges")
# for (lib in LibListx) {
#   x <- try(library(package = as.character(lib), character.only = T))
#   if (class(x) == "try-error") {
#     BiocManager::install(pkgs = lib)
#   }
# }

# load packages  --------------------------------------------------------------------
library(Rcpp)
library(devtools)
library(methods)
library(readr)
library(Biostrings)
library(plyr)
library(dplyr)
library(data.table)
library(stringr)
library(DECIPHER)
library(roperators)

# Body  --------------------------------------------------------------------
needs_loc <- "/home/neri/Downloads/hh-suite/build/bin:/home/neri/Downloads/hh-suite/build/scripts:/home/neri/bin/:/usr/local/sbin/:/usr/local/bin/:/usr/sbin:/usr/bin:/home/neri/Downloads/mmseqs/bin/:/home/neri/bbmap/"
Sys.setenv(FASTA_4_MAFFT = "/home/neri/bin/fasta36")
Sys.setenv(PATH = paste0(Sys.getenv("PATH"), needs_loc))
mod_GENETIC_CODE <- GENETIC_CODE
mod_GENETIC_CODE[which(mod_GENETIC_CODE == "*")] <- "X"
THREADS <- data.table::getDTthreads()

###### Misc. functions ######
p0 <- function(...) {
  paste0(...)
} # short call to paste0 (sep="").

len <- function(...) {
  length(...)
} # short call to length.

ClearFunct <- function() { # Remove all function from Global environment (only?).
  # Oldname: clearfunct
  list1 <- setdiff(ls(), lsf.str())
  list2 <- setdiff(ls.str(), list1)
  rm(list = list2)
}

ClearTmps <- function() { # Remove all function from Global environment (only?).
  rm(list = intersect(ls(), c("w3","w4","w1","w2","TMP2", "TMP3", "TMP4", "TMP5", "TMP_gr", "x", "testttt", "tmmp", "TMP1", "ij", "cntr", "j", "Nlb", "Nyd", "Rmems1", "RTsNodes", "tmplabs", "tmpkids", "tmpppp", "wx", "i", "w01", "w00", "gbtrees", "iz", "ix", "iy", "dfdt", "btrees", "b", "textast", "textastDF", "WNewNodTaxDF_master", "tmpclsts2", "tmpclsts", "tmppp", "tmpsubtree", "tmpdf_dist", "tmpcount", "x", "lib", "LibListx", "cx", "cx1", "xcx", "lsfiles", "faa_df", "IDFT4print", "MtmpNewNodTaxDF_toc", "tmpNewNodTaxDF_blnk", "tmpNewNodTaxDF_toc", "hits_df1", "dupdf", "aaa", "ALL_nuc.faa", "ALL_nuc.faa.dfdt", "ALL_nuc_3007.dfdt", "ALL_nuc_3007.fasta", "ALL_nuc.faa_trimmed", "workies", "Wnrws", "ToCorr", "study", "subject_was", "susrows", "Ntax", "NewTips", "mts_in_this_study", "min_len", "max_len", "www", "www1", "www2", "WWWusage", "tmptmpdf", "tmptmptmpdf", "tmptmpdf", "tmptmp3", "tmpranks", "tmplsit", "Rmems1", "Rmems2", "Rcls90", "RNC90", "Rnuc.cls9595", "wx", "w0", "w00", "w1", "ir", "lambda", "mems1", "NewNodTaxDF_template", "ir", "maxlier", "brnm", "LCARnkCntr", "AllwRank", "colstst", "AllTrees", "db2c.2", "colstst1", "distss", "exncol", "ncolneeded", "nrowsneeded", "w3", "WantedKids", "UnwantedKids", "i", "M_c", "M_t", "ntax", "l", "drows11", "drows4", "workdf_certain", "workdf1", "workDF1", "workdf2", "workDF2", "TaxMotCount", "TaxMotCount3", "TaxMotCount3.2", "RsomeDF1", "rdrp_v301.c210209.id.tab", "Parsed_RBS_motifs_gen11", "Parsed_RBS_motifs_gen4", "p3", "p4", "p5", "p6", "Orig_tmp.210112a.phylo", "fsomeDF1", "someDF", "someDF1", "tmpfasta", "tmp_4sureenaes", "w00", "w0000", "w0001", "NewNode", "tmpfilter", "testmmmm", "tmpTaxDF", "tmpmp", "tmpmat", "tmp_allnaes", "tmpa", "tmpclstmems4count", "y", "WnPrDf", "WnTxDf", "W_NNTFM", "tmperdf", "tmpdf", "tmp_notsure", "tmpNodetbl_test2", "tmprank", "tmptpmpmpmpm", "tmptreeque", "tmptreeque", "temp_file", "rdrp_v301.05.hit.tab", "RdRps3007.seqsdf", "olog", "temp2", "test13", "test23", "test1", "xcl", "tempdt", "tmpali", "tmptree", "tmptest133", "w", "AAA", "cc2", "asdsad", "adgtmpali", "mems", "mem", "cgtables", "outdf", "Nucs", "nwolf", "newtree", "rv20_v2", "ggNeoTree", "c210112a.seqsdf", "nucls", "RdRps", "bspkgs.nbu", "bspkgs.nb", "command", "kids", "Comscaf", "cpath", "reps1", "FlgType", "InrNodes2eval", "Tips2eval", "tmpvar2", "wd", "origdir", "param", "outfile", "clsts", "AsType", "Comsdf", "ya20_JAAOEH010001231_1", "omit.internal", "infile", "Taxa", "ALL_nuc.fasta", "clsts.nuc.newids.9595", "compdf", "AllIDs_vr0520", "AllIDs_set2n0_ALL_nuc", "ALL_nuc.fasta_df", "AllLCAKids", "AllIDs", "clsts.nuc.newids.9595", "AllIDs_vr0520_test", "ALL_nuc.fasta_df1", "cls9590", "info", "ModifiedW0", "rdrp_v301_162507", "list2", "states", "taxid2lineage", "taxid2lineageTrim", "cntr", "ccls", "logvar", "trimifo", "trimito", "wp", "InrNodes2Rem", "stated", "TreeQue", "maxclade", "wy", "wx", "ix", "iy", "i", "j", "cc", "a", "clade", "clades", "Yith", "izh", "label", "tmpNotReady", "tmptest1", "logvarmVizTree", "w1", "X", "x", "x1", "x2", "x3", "w0", "tree", "node", "GenDF", "cols2test", "subtaxtree", "splitedtax", "AllanoY", "Allano", "Allanotips", "Allanox", "tmplst", "tmpsum", "test1021", "test1", "test2", "test3", "test1121", "maxval", "maxX", "maxF", "LCA", "nodesdone", "nodesleft", "Nodes4Pies", "pb", "NODES2DROP", "NID", "problomatic", "thiscladelca", "tmpanc", "listerrs", "maxcladelca", "Max_Rank", "include.self", "Tips", "Tips2Rem", "cols", "w2", "w8", "mll", "nnid", "gcode", "cpath2", "cpath", "TreeQue_NoTips", "ProtectedLabs")))
  gc()
}

Odd <- function(x){
  x%%2 == 1
}

Even <- function(x){
  x%%2 != 1
}

pad <- function(str, pad = " ", width = floor(0.9 * getOption("width")), side = c("left", "right", "both"), use_length = FALSE) {
  stringi::stri_pad(str, width, side, pad, use_length)
}

TabUnroll <- function(dfdt, sep = ",", colnym = "mems", NewColNym = "values") {
  dfdt$TmpLstCol <- apply(X = dfdt, MARGIN = 1, FUN = function(x) unlist(unname(str_split(string = x[colnym], pattern = coll(sep, ignore_case = FALSE, locale = "en"), n = Inf,simplify=T))))
  dfdt3 <- with(dfdt, {
    data.frame(lapply(`dropcols<-`(dfdt, "TmpLstCol"), rep, times = lengths(TmpLstCol)), TmpCol1 = unlist(TmpLstCol))
  })
  dfdt3 <- `dropcols<-`(dfdt3, "TmpLstCol")
  dfdt3 <- Rename1Col(dfdt3, "TmpCol1", NewColNym)
  return(dfdt3)
}

Change2Date<-function(z) { # From Leah, use like this: timevec<-as.numeric(difftime(date1,date2),units="weeks")
  x=grep("date",colnames(z),ignore.case =T)
  for ( i in (1:length(x))) {
    z[ ,x[i]]=as.Date(z[ ,x[i]],"%d/%m/%Y")
  }
  z
}

QuickLoad <- function(THREADS = THREADS) {
  lsfiles <- file.info(dir("Rdatacatalog/", full.names = TRUE, pattern = "*.RDataFS"))
  fastSave::load.lbzip2(file = rownames(lsfiles)[which.max(lsfiles$mtime)], n.cores = THREADS)
}

fun_color_range <- colorRampPalette(c("#FFF2B9", "red"))   # Apply colorRampPalette From StackOverflow


QuickSave <- function(THREADS = THREADS, rcpatn = "Rdatacatalog/", name = "KvPhylo") {
  name <- "KvMeta"
  rcn <- p0("Rdatacatalog/", name, ".", Sys.Date(), ".RDataFS")
  fastSave::save.image.lbzip2(rcn, n.cores = THREADS)
}

QuickRDSSave <- function(object = p71,n.cores=THREADS,fname="tmp",returnfname=F) {
  fname <- p0(fname,".",gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2]))),".RDS")
  fastSave::saveRDS.lbzip2(n.cores = n.cores,object = gginnards::drop_vars(object),file = fname)
  print(p0("RDS saved  -  ",fname))
  if(returnfname){return(fname)}
}

QuickRDSRead <- function(n.cores=THREADS,fname="tmp") {
  fastSave::readRDS.lbzip2(n.cores = THREADS,file = fname)
}

as.char <- function(...) {
  as.character(...)
}

as.num <- function(...) {
  as.numeric(...)
}

AllTrue <- function(x) { # Are all items of a logical vector == TRUE?
  if (length(which(x == F)) == 0) {
    return(T)
  }
  return(F)
}

AllFalse <- function(x) { # Are all items of a logical vector == False?
  if (length(which(x != F)) == 0) {
    return(T)
  }
  return(F)
}

AnyFalse <- function(x) { # Are all items of a logical vector == False?
  if (length(which(x == F)) > 0) {
    return(T)
  }
  return(F)
}


show_me_the_contig <- function(full_name) { # To view on IMG/MER site.
  brkt <- (strsplit(as.character(full_name), split = "_", fixed = T))
  mt <- brkt[[1]][1]
  scf <- paste0(as.character(brkt[[1]][2]), "_", as.character(brkt[[1]][3]))
  bu_wd <- paste0("https://img.jgi.doe.gov/cgi-bin/mer/main.cgi?section=MetaScaffoldDetail&page=metaScaffoldDetail&taxon_oid=", mt, "&scaffold_oid=", scf, "&data_type=assembled")
  browseURL(bu_wd)
  invisible(bu_wd)
}

fread0 <- function(...) { # Shortened fread
  fread(quote = F, header = T, sep = "\t",na.strings = "", ...)
}


which.na <- function(x) {
  which(is.na(x))
}

which.Not.na <- function(x) {
  which(!is.na(x))
}

which.duplicated <- function(x, returnValue = F, returnalldups = T) {
  w1 <- which(duplicated(x))
  if (!returnValue) {
    if (returnalldups) {
      return(which(x %in% x[w1]))
    }
    return(w1)
  }
  if (returnValue) {
    return(unique(x[which(x %in% x[w1])]))
  }
}

# which.Not.duplicated <- function(x, returnValue = F) {
#   w1 <- which(duplicated(x))
#   if(length(w1) == 0){
#     return(seq_along(x))
#   }
#   if (!returnValue) {
#       return(!which(x %in% x[w1]))
#   }
#   if (returnValue) {
#     return(unique(x[!which(x %in% x[w1])]))
#   }
# }

whd <- function(x, ...) {
  which.duplicated(x, ...)
}
# 
# wahd <- function(x, ...) {
#   which.Not.duplicated(x, ...)
# }

wna <- function(x, ...) {
  which.na(x, ...)
}

wana <- function(x, ...) {
  which.Not.na(x, ...)
}

SharedCols <- function(df1, df2) {
  intersect(colnames(data.frame(df1)), colnames(data.frame(df2)))
}

UnSharedCols <- function(df1, df2) {
  unique(setdiff(colnames(data.frame(df2)), colnames(data.frame(df1))),setdiff(colnames(data.frame(df1)), colnames(data.frame(df2))))
}


wh <- function(x, ...) {
  which(x, ...)
}

UpdateNotation <- function(indf = IDFT, O2NDF = O2NDF, Ntype = "All") {
  if (!(Ntype %in% c("Nuc", "ORFs", "C90", "C95", "RdRp", "All"))) {
    print("Ntype most be one of: Nuc, ORFs, C90, C95, All, RdRp")
  }
  if (Ntype == "ORFs") {
    indf <- merge(Rename1Col(indf, "ORFID", "OLD_ORFID"), O2NDF, by = intersect(colnames(O2NDF), colnames(indf)), all.x = T, all.y = F)
    return(`dropcols<-`(indf, c("OLD_ORFID", "C90", "C95", "RID", "new_nuc_id", "nuc.cls9090", "cls9595_id", "RdRp_ID")))
  }
  indf <- merge(indf, O2NDF, by = intersect(colnames(O2NDF), colnames(indf)), all.x = T, all.y = F)
  if (Ntype == "Nuc") {
    return(`dropcols<-`(indf, setdiff(colnames(O2NDF), "ND")))
  }
  if (Ntype == "C90") {
    return(`dropcols<-`(indf, setdiff(colnames(O2NDF), "NC90")))
  }
  if (Ntype == "C95") {
    return(`dropcols<-`(indf, setdiff(colnames(O2NDF), "NC95")))
  }
  if (Ntype == "RdRp") {
    w1 <- which.na(indf$RID)
    if (length(w1) != 0) {
      indf$RID[w1] <- indf$RdRp_ID[w1]
    }
    return(`dropcols<-`(indf, setdiff(colnames(O2NDF), "RID")))
  }
  if (Ntype == "All") {
    return(`dropcols<-`(indf, c("nuc.cls9090", "cls9595_id", "RdRp_ID", "new_nuc_id")))
  }
}

zstdR <- function(files = files, zstdout="tmp.zstd",compression_level=22){
  system2(p0("zstd -",compression_level," -o ", zstdout," ",collapse(files)))
}



###### DF/DT functions ######
CalcPcoverage <- function(hit_df) { # Calculate profile coverage, based on alignment coordinates and profile length.
  # Oldname: calc_pCoverage
  if (AnyFalse(hit_df$p2 > hit_df$p1)) {
    tmplsit <- hit_df$p1 > hit_df$p2
    hit_df$pCoverage <- NA
    hit_df$pCoverage[tmplsit] <- (hit_df$p1[tmplsit] - hit_df$p2[tmplsit] + 1) / hit_df$pL[tmplsit]
    hit_df$pCoverage[!tmplsit] <- (hit_df$p2[!tmplsit] - hit_df$p1[!tmplsit] + 1) / hit_df$pL[!tmplsit]
    return(hit_df)
  }
  hit_df$pCoverage <- (hit_df$p2 - hit_df$p1 + 1) / hit_df$pL
  return(hit_df)
}

CalcQcoverage <- function(hit_df) { # Calculate Query coverage, based on alignment coordinates and profile length.
  if (AnyFalse(hit_df$q2 > hit_df$q1)) {
    tmplsit <- hit_df$q1 > hit_df$q2
    hit_df$qCoverage <- NA
    hit_df$qCoverage[tmplsit] <- (hit_df$q1[tmplsit] - hit_df$q2[tmplsit] + 1) / hit_df$qL[tmplsit]
    hit_df$qCoverage[!tmplsit] <- (hit_df$q2[!tmplsit] - hit_df$q1[!tmplsit] + 1) / hit_df$qL[!tmplsit]
    return(hit_df)
  }
  hit_df$qCoverage <- (hit_df$q2 - hit_df$q1 + 1) / hit_df$qL
  return(hit_df)
}

`dropcols<-` <- function(df, dcols) { # Replacement function to remove columns.
  df[, dcols] <- NULL
  df
}
`droprows<-` <- function(df, drows) { # Replacement function to remove columns.
  rwnms <- rownames(df)
  if (!is.numeric(drows)) {
    w <- which(rwnms %in% drows)
    return(`droprows<-`(df, w))
  }
  if (length(drows) > 0) {
    df <- df[-drows, ]
    rownames(df) <- rwnms[-drows]
    return(df)
  }
  df
}

CompareDFDT <- function(dfdt1, dfdt2) {
  # base '==' doesn't work for DFs with nested lists?
  if (nrow(dfdt1) != nrow(dfdt2) || ncol(dfdt1) != ncol(dfdt2)) {
    return(F)
  }
  logdf <- data.frame(nrow = nrow(dfdt1), ncol = ncol(dfdt2))
  for (Ir in 1:nrow(dfdt1)) {
    for (Ic in 1:ncol(dfdt2)) {
      logdf[Ir, Ic] <- (unlist(dfdt1[Ir, Ic]) == unlist(dfdt2[Ir, Ic]))
    }
  }
  logvar <- apply(logdf, MARGIN = 2, FUN = function(x) unique(x))
  if (length(which(unlist(logvar) == F)) > 0) {
    return(F)
  }
  return(T)
}

RemoveEmptyStrRow <- function(dfdt) {
  # Oldname: rename1col
  w1 <- unlist(apply(dfdt, MARGIN = 1, FUN = function(x) AllTrue(x == "")))
  if (length(w1) != 0) {
    if (length(which(w1)) != 0) {
      dfdt <- `droprows<-`(dfdt, drows = which(w1))
    }
  }
  dfdt
}

RemoveEmptyStrCol <- function(dfdt) {
  w1 <- unlist(apply(dfdt, MARGIN = 2, FUN = function(x) AllTrue(x == "")))
  if (length(w1) != 0) {
    if (length(which(w1)) != 0) {
      dfdt <- `dropcols<-`(dfdt, dcols = which(w1))
    }
  }
  dfdt
}

RemoveNACols <- function(dfdt) {
  w1 <- unlist(apply(dfdt, MARGIN = 2, FUN = function(x) AllTrue(is.na(x))))
  if (length(w1) != 0) {
    if (length(which(w1)) != 0) {
      dfdt <- `dropcols<-`(dfdt, dcols = which(w1))
    }
  }
  dfdt
}

RemoveNARows <- function(dfdt){
  w1 <- unlist(apply(dfdt, MARGIN = 1, FUN = function(x) AllTrue(is.na(x))))
  if (length(w1) != 0) {
    if (length(which(w1)) != 0) {
      dfdt <- `droprows<-`(dfdt,which(w1))
    }
  }
  dfdt
}


Rename1Col <- function(dfdt, colnm, newcolnm) {
  # Oldname: rename1col
  if (is.matrix(dfdt)) {
    dfdt <- data.frame(dfdt, stringsAsFactors = F)
  }
  if (is.numeric(colnm)) {
    names(dfdt)[colnm] <- newcolnm
    return(dfdt)
  }
  names(dfdt)[names(dfdt) == colnm] <- newcolnm
  dfdt
}
RemovesSparseols <- function(indf) {
  # Oldname: removesparsecols
  indf <- data.frame(indf)
  rembols <- apply(indf, 2, unique)
  nurembols <- simplify2array(lapply(rembols, FUN = length))
  return(indf[, which(nurembols != 1)])
}
Trimifo <- function(x, y) {
  # Oldname: trimifo
  return(max(1, as.numeric(x["q1"]) - (y * (as.numeric(x["p1"]) - 1))))
}
Trimito <- function(x, y) {
  # Oldname: trimito
  return(min(as.numeric(x["qL"]), as.numeric(x["q2"]) + (y * (as.numeric(x["pL"]) - as.numeric(x["p2"])))))
}
HeaderBreaker <- function(input_df, sep = ".", clb = "id") {
  # Oldname: Header_breaker_gen
  breaked <- within(input_df, splitted = data.frame(do.call("rbind", strsplit(as.character(as.data.frame(input_df)[, c(clb)]), sep, fixed = T)), stringsAsFactors = F))
  return(breaked)
}
HeaderBreakerCb <- function(input_df, sep = ".", clb = "id", nclb = hedclbs) {
  i <- colnames(input_df)
  input_df <- cbind(input_df, data.frame(do.call("rbind", strsplit(as.character(as.data.frame(input_df)[, c(clb)]), sep, fixed = T)), stringsAsFactors = F))
  colnames(input_df) <- c(i, nclb)
  return(input_df)
}
WriteWolfTbl <- function(Rtbl, filepath, ...) {
  write.table(Rtbl, filepath, quote = F, row.names = F, col.names = T, sep = "\t", ...)
}

fwrite0 <- function(sep='\t',file,x,quote = T,verbose = T,...){
  fwrite(sep = sep, file = file, x = x,quote = quote,verbose = verbose,...)
}

WriteXlsx <- function(dfdt, filepath = "tmpout.xlsx", ...) {
  writexl::write_xlsx(x = dfdt, path = filepath, col_names = T, format_headers = T, ...) #
}

ReadXlsx <- function(filepath = "tmpout.xlsx", ...) {
  distinct(readxl::read_xlsx(filepath, na = "NA", col_types = "text")) #
}

MergeXtYf <- function(Xdf, Ydf, ...) {
  merge(x = Xdf, y = Ydf, all.x = T, all.y = F, ...)
}

###### BioFormats / Genomic Ranges functions ######
readGFF <- function(...) {
  rtracklayer::readGFF(...)
}

CreateCSVFromGBFF <- function(input.gbff) { # Forked from https://github.com/federiva/Monkeys-Working
  # The .sh script must be saved in the working directory
  path.to.script <- "/home/neri/Documents/GitHub/Monkeys-Working/mitochondrial_DNA/r_project_mitochondrial_DNA/bash_scripts/CreateCSVFromGBFF.sh"
  arguments.to.bash <- paste(path.to.script, input.gbff, sep = " ")

  system2(command = "/bin/bash", args = arguments.to.bash)
  print("The resulting .csv file is saved in the working
        directory as indexes_summary.csv")
}

GetNCBITaxonomyChildNodes <- function(nodes.dmp.path, target.taxid) { # Forked from https://github.com/federiva/Monkeys-Working
  # nodes.dmp.path is a STRING which is the path to the nodes.dmp
  #  file downloaded from NCBI:
  #  ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdmp.zip
  # target.taxid is a  STRING which is the upper NCBI taxonomy
  #  taxid from which we want to get the child's taxids
  # The output of the function is a .csv file ("|" delimited)
  #  saved in the current working directory named child_nodes.csv
  # The function assumes that the .sh script
  #  GetNCBITaxonomyChildNodes.sh is located in the
  #  working directory.
  path.to.script <- "/home/neri/Documents/GitHub/Monkeys-Working/mitochondrial_DNA/r_project_mitochondrial_DNA/bash_scripts/GetNCBITaxonomyChildNodes.sh"
  arguments.to.bash <- paste(path.to.script, nodes.dmp.path, target.taxid, sep = " ")
  system2(command = "/bin/bash", args = arguments.to.bash)
  cat("The resulting .csv file is saved in the working directory as child_nodes.csv")
}

ExtractFastaSeqsFromGBFF <- function(path.to.folder.with.gbff.files) { # Forked from https://github.com/federiva/Monkeys-Working
  # The .sh script must be located in the working directory
  path.to.script <- "/home/neri/Documents/GitHub/Monkeys-Working/mitochondrial_DNA/r_project_mitochondrial_DNA/bash_scripts/ExtractFastaSeqsFromGBFF.sh"
  arguments.to.bash <- paste(path.to.script, path.to.folder.with.gbff.files, sep = " ")
  system2(command = "/bin/bash", args = arguments.to.bash)
  cat("The individual .fasta sequences were saved in ./gbff_filtered_files/fasta_seqs")
}


Grange2gbRecord <- function(Grng, sequences, contig_name = "Tuliv", contig_length = 24864,
                            source_name = "Tuliv", circ = "circular", DRNA = "DNA", organism = "Tuliv",
                            keywords = "Virus", datee = "10-JAN-2021", type = "PHG",
                            workdir = "/home/neri/scratch/X2gbk/",
                            outputF = "/home/neri/scratch/X2gbk/testing.gbk",
                            inputF = "/media/HDD1/uri/projects/Isra/Tuli/New/TuliV.fasta") {

  #### Fake the header (TODO:change to sprintf)
  Fheader <- p0("LOCUS       ", contig_name, "              ", contig_length, " bp    ", DRNA, "     ", circ, " ", type, " ", datee, "
DEFINITION  ", source_name, "
ACCESSION   ", contig_name, "
VERSION     ", contig_name, ".1
KEYWORDS    ", keywords, "
SOURCE      ", source_name, "
  ORGANISM  ", source_name, "
REFERENCE   1  (bases 1 to ", contig_length, ")
  AUTHORS   Lorem,I., Dolor,S., A,C.E., Sed,D.E., Tempor,I.
  TITLE     Lorem ipsum dolor sit amet, consectetur adipiscing elit
            sed do eiusmod tempor incididunt ut labore et dolore
  JOURNAL   Unpublished
FEATURES             Location/Qualifiers
     source          1..", contig_length, "
                     /organism=\"", source_name, "\"
                     /mol_type=\"genomic ", DRNA, "\" ")
  write(Fheader, outputF, )

  
  
  Grng@elementMetadata@listData[["type"]] =  list(rep("CDS",len(Grng)))
  Grng$Parent <- Grng$ID
  #### Write the features (TODO:change to sprintf)
  for (i in 1:length(Grng@ranges)) {
    print(i)
    tmptype <- Grng@elementMetadata@listData[["type"]][i]
    tmpstrt <- Grng@ranges@start[i]
    tmpend <- (Grng@ranges@start[i] + Grng@ranges@width[i])
    print((tmpend - tmpstrt) %% 3)
    if (as.character(Grng@strand)[i] == "-") {
      tmpstrand <- -1L
      loci <- p0("complement(", tmpstrt, "..", tmpend - 1, ")")
    } else {
      tmpstrand <- +1L
      loci <- p0(tmpstrt, "..", tmpend)
    }
    if (tmptype == "gene") {
      cat((p0("     ", tmptype, "            ", loci)), file = outputF, sep = "\n", append = T)
      cat(p0(gsub(toString(rep(" ", 11)), pattern = ",", replacement = ""), "/locus_tag=\"", Grng[i]$ID, "\""), file = outputF, sep = "\n", append = T)
    }
    if (tmptype == "CDS") {
      cat((p0("     ", tmptype, "             ", loci)), file = outputF, sep = "\n", append = T)
      cat(p0(gsub(toString(rep(" ", 11)), pattern = ",", replacement = ""), "/locus_tag=\"", Grng[i]$Parent, "\""), file = outputF, sep = "\n", append = T)
      cat(p0(gsub(toString(rep(" ", 11)), pattern = ",", replacement = ""), "/product=\"", Grng[i]$Name, "\""), file = outputF, sep = "\n", append = T)
      cat(p0(gsub(toString(rep(" ", 11)), pattern = ",", replacement = ""), "/protein_id=\"", Grng[i]$ID, "\""), file = outputF, sep = "\n", append = T)
    }
  }
  #### Write the sequence (Adapted from biofiles).
  if (length(seq = sequences) > 0L) {
    lineno <- seq(from = 1, to = seq@ranges@width, by = 60)
    lines <- seq_along(lineno)
    n_lines <- length(lines)
    s <- character(n_lines)
    for (i in lines) {
      seqw <- ifelse(i < n_lines, i * 60, seq@ranges@width)
      seqs <- XVector::toString(XVector::subseq(seq, 1 + (i - 1) * 60, seqw))
      nnn <- seq(1, nnncc = nchar(seqs), by = 10)
      s[i] <- paste0(substring(seqs, nnn, c(nnn[-1] - 1, nnncc)), collapse = " ")
    }
    s <- sprintf("%+9s %s", lineno, s)
    cat("\nORIGIN", file = outputF, sep = "\n", append = T)
    cat(s, file = outputF, sep = "\n", append = T)
    cat("//", file = outputF, append = T)
  } else {
    cat("\n//", file = outputF, append = T)
  }
}

###### Tree functions ######
Tree2Edgetbl <- function(tree) {
  # Consider replacnig with tidytree::as_tibble(tree)
  EdgeTbl <- data.frame(
    "parent_node" = tree$edge[, 1],
    "node" = tree$edge[, 2])
  EdgeTbl$node.label = unlist(NID2NLabel(tree,EdgeTbl$node))
  return(EdgeTbl)
}

QuickTreePDF <- function(name="QuickTree.pdf"){
  Name <- gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2])))
  pdf(
    file = p0("TmpTree.",Name,".",name),
    width = 50,
    height = 50,
    title = "ArgMax with CRISPR and Lysis matches tips",
    useDingbats = T,
  )
}

QuickTreePDF2 <- function(name="QuickTree.pdf"){
  Name <- gsub(pattern = ", ",replacement = ".",x = toString.default(gsub(pattern = ":",replacement = "-",x = str_split_fixed(string =  Sys.time(),pattern = " ",n = 3)[,1:2])))
  cairo_pdf(antialias = "none",
    filename = p0("TmpTree.",Name,".",name),
    width = 50,
    height = 50)#,
    # title = "ArgMax with CRISPR and Lysis matches tips",
  # )
}

RepositionClade <- function(tree, clade, where, position = NULL) {
  # drops the clade from the tree, then binds it in where and position
  cladetree <- extract.clade(phy = tree, node = clade, collapse.singles = F)
  tree <- treeio::drop.tip(object = tree, tip = phangorn::Descendants(tree, clade, type = "tips")[[1]])
  tree <- bind.tree(tree, cladetree, where = where, position = position)
  return(tree)
}

isTip <- function(tree, node, ...) {
  if(class(node) == "character"){
    return(node %in% tree$tip.label)
  }
  if(class(node) == "numeric"){
    return(node <= Ntip(tree))
  }
    # treeio::isTip(.data = tree, .node = node, ...)
}

#### Keep forgetting which package has this functions, so binding them to more memrobale function names.
GetKids <- function(tree, node, Return_Labels = F, ...) {
  if (!Return_Labels) {
    return(phangorn::Descendants(x = tree, node = node, ...))
  }
  return(lapply(phangorn::Descendants(x = tree, node = node, ...), function(x) NID2NLabel(tree, x)))
}
GetParents <- function(tree, node, ...) {
  phangorn::Ancestors(x = tree, node = node, ...)
}

DropNodeKids <- function(tree = WorkTree, node) { # 2022 Version.2
  if(class(node) == "character"){
    Nodelabel <- node
    nodeid <- NLabel2NID(tree,node)
  }
  if(class(node) == "numeric"){
    nodeid <- node
    Nodelabel <- NID2NLabel(tree,node)
  }
  w1 <- wh(tree$edge[,2] == nodeid)
  Len2Parent <- tree$edge.length[w1]
  ParentLabel <- NID2NLabel(tree,tree$edge[w1,1])
  Tips2Rem <- unlist(GetKids(tree, node = nodeid, type = "all",Return_Labels =T))
  tree <- ape::drop.tip(phy = tree, tip = Tips2Rem, trim.internal = T,collapse.singles = F)
  tree <- phangorn::add.tips(tree = tree, tips = Nodelabel, where = NLabel2NID(tree,ParentLabel), edge.length = Len2Parent)
  return(tree)
}

GetUnwantedKids <- function(tree, node, omit.internal = T, asids = F, ...) {
  node <- as.integer(na.omit(node))
  Tips2eval <- tree$tip.label[unique(unlist(GetKids(tree, node, type = "tips")))]
  InrNodes2eval <- tree$node.label[setdiff(unique(unlist(GetKids(tree, node, type = "all"))), 1:Ntip(tree)) - Ntip(tree)]
  # ProtectedLabs=tree$node.label[node-Ntip(tree)]
  # if (include.self==F){
  #   InrNodes2eval=setdiff(InrNodes2eval,ProtectedLabs)
  # }
  Tips2eval <- unique(union(Tips2eval, InrNodes2eval))
  LCA <- castor::get_mrca_of_set(tree = tree, descendants = Tips2eval)
  if (omit.internal) {
    AllLCAKids <- tidytree::nodelab(tree = tree, GetKids(tree, LCA, type = "tips")[[1]])
  }
  if (!omit.internal) {
    AllLCAKids <- tidytree::nodelab(tree = tree, GetKids(tree, LCA, type = "all")[[1]])
  }
  UnwantedKids <- setdiff(AllLCAKids, Tips2eval)
  if (asids) {
    UnwantedKids <- NLabel2NID(tree, UnwantedKids)
  }
  return(UnwantedKids)
}

SplitTree2MonoPhyloByTips <- function(tree, node) {
  # EDIT DESCRIPTION!
  # TreeTools:: input tree searched for LCA of given nodes/tips, then, if the input tips define a monophyletic group returns a subtreee
  # with that LCA as root, ommitting all tips and nodes not in the monophyletic group. If the group of tips doesn't define a monophyletic group, split it into N subgroups that do define monophyletic groups with those tips.
  # Note! tree needs to have labels (tips and nodes) because I can't be bothered right now.
  LCA <- get_mrca_of_set(tree, node)
  tmpsubtree <- get_subtree_at_node(tree = tree, node = (LCA - Ntip(tree)))$subtree
  NewNode <- NLabel2NID(tmpsubtree, NID2NLabel(tree, node))

  UnwantedKids <- unlist(GetUnwantedKids(tree = tmpsubtree, node = NewNode, asids = T))
  WantedKids <- unique(unlist(GetKids(tree = tmpsubtree, node = NewNode, type = "tips")))
  tmptreeque <- get_tree_traversal_root_to_tips(tmpsubtree, include_tips = T)$queue

  tmpdf <- data.frame("Node" = c(1:(tmpsubtree$Nnode + Ntip(tmpsubtree))), "Node_Label" = c(tmpsubtree$tip.label, tmpsubtree$node.label), "Kids" = NA, "Has_unwanted" = T, "muster" = F, "IsTip" = F)
  tmpdf$Kids <- GetKids(tree = tmpsubtree, node = tmpdf$Node, type = "all")
  tmpdf$Nunwanted <- unlist(mclapply(tmpdf$Node, FUN = function(x) unlist(length(setdiff(unlist(tmpdf[x, "Kids"]), unlist(WantedKids)))), mc.preschedule = T, mc.set.seed = T, mc.silent = F, mc.cores = THREADS, mc.cleanup = F, mc.allow.recursive = T, affinity.list = NULL))
  tmpdf$Has_unwanted <- unlist(mclapply(tmpdf$Node, FUN = function(x) {
    (if (tmpdf[x, "Nunwanted"] == 0) {
      F
    } else {
      T
    })
  }, mc.preschedule = T, mc.set.seed = T, mc.silent = F, mc.cores = THREADS, mc.cleanup = F, mc.allow.recursive = T, affinity.list = NULL))
  tmpdf$Nwanted <- unlist(mclapply(tmpdf$Node, FUN = function(x) (length(intersect(unlist(tmpdf[x, "Kids"]), unlist(WantedKids)))), mc.preschedule = T, mc.set.seed = T, mc.silent = F, mc.cores = THREADS, mc.cleanup = F, mc.allow.recursive = T, affinity.list = NULL))
  tmpdf$Has_wanted <- unlist(mclapply(tmpdf$Node, FUN = function(x) {
    (if (tmpdf[x, "Nwanted"] == 0) {
      F
    } else {
      T
    })
  }, mc.preschedule = T, mc.set.seed = T, mc.silent = F, mc.cores = THREADS, mc.cleanup = F, mc.allow.recursive = T, affinity.list = NULL))
  tmpdf <- filter(tmpdf, Has_wanted == T)
  tmpdf <- filter(tmpdf, Has_unwanted == F)

  tmpdf$IsTip[which(tmpdf$Node <= Ntip(tmpsubtree))] <- T
  kidsleft <- as.integer(WantedKids)
  kidsdone <- c()
  tmptreeque_rel <- intersect(tmptreeque, tmpdf$Node)
  for (w in (tmptreeque_rel)) {
    print(length(kidsleft))
    # print(w)
    j <- which(tmpdf$Node == w)
    # if(tmpdf$Has_unwanted[j]){
    #   next
    # }
    leny <- length(intersect(unlist(tmpdf$Kids[j]), kidsleft))
    if (leny > 0) {
      kidsdone <- unique(union(kidsdone, as.integer(unlist(tmpdf$Kids[j]))))
      tmptreeque_rel <- setdiff(tmptreeque_rel, kidsdone)
      kidsleft <- setdiff(kidsleft, kidsdone)
      tmpdf$muster[j] <- T
    }
    if (leny == 0) {
      tmpdf$muster[j] <- F
    }
    if (length(kidsleft) == 0) {
      break
    }
  }
  return(filter(tmpdf, muster))
}

NLabel2NID <- function(tree, label) {
  tidytree::nodeid(tree, label)
}

NID2NLabel <- function(tree, tip) {
  if(len(tip) == 1){
    if(tip <= Ntip(tree)){
      return(tree$tip.label[tip])
    }
    return(tree$node.label[tip - Ntip(tree)])
  }
  if(len(tip) > 1){
    sapply(X = tip, FUN = function(x) NID2NLabel(tree,x))
  }
}


###### B/X string functions ######
XString2DF <- function(faa, input_was = "new_name", trimwhite = F, seqcolname = "seq", addlength = T) {
  # Oldname: XString2df
  o1 <- as.character(faa)
  outdf <- data.frame(xyz = names(o1), "seq" = unname(o1), Length = faa@ranges@width)
  if (trimwhite == T) {
    outdf$xyz <- trimws(outdf$xyz)
  }

  colnames(outdf) <- c(toString(input_was), seqcolname, "Length")
  if (addlength == F) {
    outdf$Length <- NULL
  }
  return(outdf)
}

DF2XString <- function(dfdt, seqcol = "seq", nmcol = "id") {
  # Oldname: df2XString
  dfdt <- as.data.frame(dfdt)
  outxstring <- BStringSet(x = dfdt[, seqcol])
  names(outxstring) <- dfdt[, nmcol]
  return(outxstring)
}

XString2KmerDT <- function(infaa = BS_all_mot_sssp, abcde = c("H", "E", "C"), k = 3, out_probablity = F, input_was = "name_mot") { # out_probablity == return output as absoulte number (count) or as probablity
  allKmers <- BStringSet(mkAllStrings(abcde, k, fast.moving.side = "right"))
  names(allKmers) <- allKmers
  Kmer_counts <- data.table(t(vcountPDict(allKmers, infaa[], max.mismatch = 0, min.mismatch = 0, with.indels = F, fixed = T, algorithm = "auto", verbose = T)))
  Kmer_counts <- cbind(c(names(infaa)), Kmer_counts)
  colnames(Kmer_counts) <- c(input_was, names(allKmers))
  Kmer_counts <- Kmer_counts[, lapply(.SD, as.numeric), by = input_was]
  mertypes <- length(allKmers)
  if (out_probablity == F) {
    return(Kmer_counts)
  } else {
    # Kmer_counts_1=Kmer_counts
    Kmer_counts$sumr <- as.numeric(apply(Kmer_counts[, -1], 1, FUN = function(x) sum(x)))
    for (mer in names(allKmers)) { # Figure out a better way to do this
      Kmer_counts[, mer] <- apply(Kmer_counts, 1, FUN = function(x) (as.numeric(x[mer])) / as.numeric(x["sumr"]))
    }
    return(`dropcols<-`(Kmer_counts, "sumr"))
  }
  # Need to refactor in O(n*(log(n))), get only exisiting Kmers, string by string, added to Dict every kmer, +1 to the count at string name X that kmer)
}
AAcoor2NAcoor <- function(frm, AAc1, AAc2, qL) { # Transform alignment coordinates from protein locations to ~locationn on the (coding?) nucleic acid sequence.
  # print(frm)
  frm  <- as.numeric(unlist(unname(trimws(frm))))
  AAc1 <- as.numeric(unlist(unname(trimws(AAc1))))
  AAc2 <- as.numeric(unlist(unname(trimws(AAc2))))
  qL   <- as.numeric(unlist(unname(trimws(qL))))
  
  if (frm > 0) {
    NAc1 = frm + AAc1*3 - 3
    NAc2 = frm + AAc2*3 - 1
  }
  if (frm < 0) {
    NAc1 = qL + frm - AAc1*3 + 4
    NAc2 = qL + frm - AAc2*3 + 2
  }
  
  return(p0(NAc1,"..", NAc2))
}

Prodigal_AA2NA <- function(ORFlen, ORFstart_nuc, ORFend_nuc, AA1, AA2,Strand) { # Transform alignment coordinates from a prodgial or similar predicted region. Not suitable for 6frxs!.
  AA1 <- as.num(AA1)
  AA2 <- as.num(AA2)
  ORFlen <- as.num(ORFlen) 
  ORFstart_nuc <- as.num(ORFstart_nuc)
  ORFend_nuc <- as.num(ORFend_nuc)
  if(Strand == "-"){
    DeltaQ <- ORFlen + 1 - (AA2 + 1)
    NA2 <- (DeltaQ * 3) + ORFstart_nuc
    NA1 <- NA2 + (3 * (AA2 - AA1 -1))
  }
  if(Strand == "+"){
    NA1 <- (((AA1) - 1) * 3) + (ORFstart_nuc)
    NA2 <- ORFstart_nuc + ((AA2 - 1) * 3)
  }
  return(p0(NA1,"..",NA2))
}

AAcoor2NAcoor_df <- function(dfdt, frmcol, AAc1col, AAc2col, outdf = T) { # Transform alignment coordinates from protein locations to ~locationn on the (coding?) nucleic acid sequence.
  # # Reverse:
  # n1 <- (q2 * 3) - frame
  # n2 <- ((q1 - 1) * 3) - frame -1
  # 
  # # Forward:
  # n1 <- ((q1 - 1) * 3) + frame
  # n2 <- (q2 * 3) + frame +1
  # 
  # width = n2 - n1 + 1
  data.table::setDT(dfdt)
  rv <- which(dfdt[, ..frmcol] < 0)
  fr <- which(dfdt[, ..frmcol] > 0)
  dfdt$NAcoor <- ""
  dfdt$NAc1 <- -1
  dfdt$NAc2 <- -1
  dfdt$NAc1[rv] <- as.integer(unlist(unname(dfdt[rv, ..AAc2col] * 3))) - as.integer(unlist(unname(dfdt[rv, ..frmcol])))
  dfdt$NAc2[rv] <- as.integer(unlist(unname((dfdt[rv, ..AAc1col] - 1) * 3))) - as.integer(unlist(unname(dfdt[rv, ..frmcol]))) -1
  dfdt$NAc1[fr] <- as.integer(unlist(unname((dfdt[fr, ..AAc1col] - 1) * 3))) + as.integer(unlist(unname(dfdt[fr, ..frmcol])))
  dfdt$NAc2[fr] <- as.integer(unlist(unname(dfdt[fr, ..AAc2col] * 3))) + as.integer(unlist(unname(dfdt[fr, ..frmcol]))) +1
  # dfdt$NAc2[which(dfdt$NAc2 > dfdt$Length)] <- dfdt$Length[which(dfdt$NAc2 > dfdt$Length)]
  # dfdt$NAc1[which(dfdt$NAc1 > dfdt$Length)] <- dfdt$Length[which(dfdt$NAc1 > dfdt$Length)]
  # dfdt$NAc1[which(dfdt$NAc1<0)]=dfdt$Length[which(dfdt$NAc1>dfdt$Length)]

  dfdt[, "NAcoor"] <- p0(dfdt$NAc1, "-", dfdt$NAc2)
  if (outdf == T) {
    return(dfdt)
  }
  return(dfdt$NAcoor)
}

# dfdt$NAc2[which(dfdt$NAc2 > dfdt$Length)] <- dfdt$Length[which(dfdt$NAc2 > dfdt$Length)]
# dfdt$NAc1[which(dfdt$NAc1 > dfdt$Length)] <- dfdt$Length[which(dfdt$NAc1 > dfdt$Length)]

AAcoor2NAcoor_dFAST <- function(dfdtrw) { # same as AAcoor2NAcoor_df but quicker and assumes specific colnames.
  dfdtrw <- data.frame(dfdtrw)
  dfdtrw$frame <- as.numeric(trimws(dfdtrw[, "frame"]))
  # AAc1=as.numeric(dfdtrw[,"q1"])
  # AAc2=as.numeric(dfdtrw[,"q2"])
  dfdtrw[,c("NAcoor")] <- apply(dfdtrw, MARGIN = 1, FUN = function(x) (AAcoor2NAcoor(frm = x["frame"], AAc1 = x["Oq1"], AAc2 = x["Oq2"],qL = x["QL"])))
  return(dfdtrw)
}

XSDNARedunFilter <- function(xsnuc) { # removes exact duplicates, with perfect matches. Returns the filtered set in no particular order (discards one of the duplicates)
  # Oldname: scaffold_redundany_filter_revcom
  xsnuc_rev <- reverseComplement(xsnuc)
  rev_reg_xsnuc <- unique(union(xsnuc, xsnuc_rev))
  uniquers <- unique(names(rev_reg_xsnuc))
  return(xsnuc[uniquers])
}

Trim2Core2 <- function(hits_df = Motif_df, faa = FL_longest_members, add_genetic_code = F, Out_faa = T, out_df = F, input_was = "new_name", subject_was = "profile", Expand_projectionX = 0) {
  # Oldname: Trim2motif_V3
  faa_df <- XString2DF(faa, input_was = input_was, seqcolname = "seq")
  workies <- intersect(hits_df[, input_was], names(faa))
  hits_df <- hits_df[which(hits_df[, input_was] %in% workies), ]
  work_faa <- faa[hits_df[, input_was]] # This'll keep duplicates AND the order.
  hits_df <- merge(hits_df, faa_df[, c(input_was, "seq")], by = c(input_was), all.x = T, all.y = F)
  # which.na(hits_df$seq)
  hits_df$QnS <- p0(hits_df[, input_was], ".", hits_df[, subject_was])
  if (Expand_projectionX != 0) {
    hits_df$extract_from <- apply(hits_df, 1, trimifo, Expand_projectionX)
    hits_df$extract_to <- apply(hits_df, 1, trimito, Expand_projectionX)
    out_faa <- BStringSet(hits_df$seq, start = hits_df$extract_from, end = hits_df$extract_to)
  } else {
    out_faa <- BStringSet(hits_df$seq, start = hits_df$q1, end = hits_df$q2)
  }
  names(out_faa) <- hits_df$QnS
  if (add_genetic_code == T) {
    hits_df$Genetic_Code <- "Standard"
    hits_df$Genetic_Code[unique(grep("X", x = out_faa[]))] <- "Non-Standard"
  }
  output <- list()
  if (Out_faa == T) {
    output <- append(output, out_faa)
  }
  if (out_df == T) {
    hits_df <- data.frame(hits_df, stringsAsFactors = F)
    output <- list(output, hits_df[, setdiff(colnames(hits_df), c("seq"))])
  }
  return(output)
}

###### Parsers and Runners ######
MMseq2_outfmt6_cols <- c("subject_name", "evalue", "gapopen", "pident", "nident", "q1", "q2", "qL", "p1", "p2", "pL", "ali_len", "raw", "score", "frame", "mismatch", "qcov", "tcov")
hmmsearh_cols <- c("r1", "qL", "subject_name", "r2", "pL", "evalue", "score", "bias", "#", "of", "c-Evalue", "i-Evalue", "score2", "bias", "p1", "p2", "q1", "q2", "env_from", "env_to", "acc", "r3", "r4")
psiblast_cols <- c("pident", "q1", "q2", "p1", "p2", "qL", "pL", "ali_len", "evalue", "score", "subject_name")
hhsearch_cols <- c("subject_name", "pCoverage", "ali_len", "pL", "mismatch", "gapOpen", "q1", "q2", "p1", "p2", "Probab", "evalue", "score")
DaimondP_cols <- c("subject_name", "evalue", "gapopen", "pident", "ali_len", "p1", "p2", "q1", "q2", "pL", "qL", "mismatch", "score")
DaimondP_rev <- c("subject_name", "pident", "ali_len", "q1", "q2", "p1", "p2", "qL", "pL", "mismatch", "gapopen", "evalue", "score")
Blastn_cols <- c("query_name", "subject_name", "pident", "qL", "pL", "q1", "q2", "p1", "p2", "ali_len", "evalue", "score")
unicols <- c("profile, qL, pL, p1, p2, q1, q2, score, evalue, ali_len, pCoverage")
MMseq2_Long_fmt6_cols <- c("query_name", "subject_name", "evalue", "gapopen", "pident", "nident", "q1", "q2", "qL", "p1", "p2", "pL", "ali_len", "raw", "score", "qframe", "mismatch", "qcov", "tcov")

colnyms <- data.frame("psiblast" = toString(psiblast_cols), "mmseqs" = toString(MMseq2_outfmt6_cols), "hhsearch" = toString(hhsearch_cols), "hmmsearch" = toString(hmmsearh_cols), "diamondp" = toString(DaimondP_cols), "uni" = unicols, stringsAsFactors = F)
GenericHitsParser <- function(inpt = "hit.tsv", versiza = versiza, Query2Profile = T, breakhdrs = T, input_was = "Scaf_ID_frame", Cull_hits = F, search_tool = "psiblast", reducecols = T, calc_pcoverage = T, colsnms = "provide", CullCol = "score", hedsep = ".", hedclbs = c("id", "frame")) {
  # Oldname: Generic_Hits_prsr_6frx
  if (!(search_tool %in% c("psiblast", "mmseqs", "hmmsearch", "hhsearch", "diamondp"))) {
    print("If search tool isn't psiblast, mmseqs, diamondp, hhmsearch or hmmsearch, provide the raw input table colnames in the arg colsnms")
    # return()
  } else {
    colsnms <- c(input_was, stringr::str_split(colnyms[1, which(colnames(colnyms) == search_tool)], ", ", simplify = T))
  }
  hits_tbl <- try(fread(inpt, stringsAsFactors = F, col.names = colsnms, sep = "\t"))
  if ((search_tool == "mmseqs")) {
    hits_tbl <- subset(hits_tbl, select = -c(6, 4, 16))
  }
  if ((search_tool == "hmmsearch")) {
    hits_tbl <- subset(hits_tbl, select = -c(2, 5, 23, 24))
    hits_tbl$ali_len <- hits_tbl$q2 - hits_tbl$q1
  }
  if (Cull_hits == T) {
    dt <- data.table(hits_tbl)
    min_dt <- data.frame(dt[, max(score), by = input_was])
    colnames(min_dt) <- c(input_was, CullCol)
    hits_tbl <- merge(min_dt, hits_tbl, by = c(input_was, CullCol), all.x = T, all.y = F, stringsAsFactors = F) # Hard culling
  }
  hits_tbl$profile <- hits_tbl$subject_name
  if (breakhdrs == T) {
    hits_tbl <- HeaderBreakerCb(input_df = hits_tbl, sep = hedsep, clb = input_was, nclb = hedclbs)
    if (reducecols == T) {
      hits_tbl <- subset(hits_tbl, select = c(input_was, hedclbs, intersect(str_split(colnyms[1, "uni"], ", ")[[1]], colnames(hits_tbl))), stringsAsFactors = F)
    }
  } else {
    if (reducecols == T) {
      hits_tbl <- subset(hits_tbl, select = c(input_was, intersect(str_split(colnyms[1, "uni"], ", ")[[1]], colnames(hits_tbl))), stringsAsFactors = F)
    }
  }
  if (calc_pcoverage == T) {
    hits_tbl <- CalcPcoverage(hits_tbl)
  }
  gc()
  return(hits_tbl)
}
ReadMcl <- function(mclfile, col1prefix = "motif") {
  # Oldname: readmcl
  cx <- scan(mclfile, what = "", sep = "\n", )
  xcx <- strsplit(cx, "[[:space:]]+") # Separate elements by one or more whitespaces...
  cluster_df <- data.frame("cls" = paste0(col1prefix, ".", 1:length(xcx)), "reps" = "", stringsAsFactors = F)
  for (i in 1:nrow(cluster_df)) {
    cluster_df$reps[i] <- (xcx[i])
  }
  return(cluster_df)
}

ReadABC <- function(ABCfile, clsnms = c("reprs", "mems")) {
  cx <- fread(input = ABCfile, data.table = F, skip = 0, header = F, sep = "\n")
  cx <- distinct(cx)
  cx <- as.data.table(stringi::stri_split_fixed(as.character(cx[, "reprs"]), pattern = " ", n = 2, simplify = T))
  cx[, lmems := list(stringi::stri_split_fixed(as.character(.SD[]), pattern = " ")), by = "V1", .SDcols = "V2"]
  colnames(cx) <- c(clsnms, "lmems")
  return(cx)
}

ReadABC2 <- function(ABCfile, clsnms = c("reprs", "mems")) {
  cx <- fread(input = ABCfile, data.table = F, skip = 0, header = F, sep = "\t")
  cx <- distinct(cx)
  xcx <- strsplit(cx$V2, "[[:space:]]+") # Separate elements by one or more whitespaces...
  cx <- as.data.table(cx)[, lmems := xcx]
  colnames(cx) <- c(clsnms, "lmems")
  return(cx)
}

ReadDBN <- function(DBNfile, includes_MFE = T, return_type = "DF") { # Read DBN file
  cx <- scan(DBNfile, what = "", sep = "\n")
  Ncx <- gsub(pattern = ">", replacement = "", x = cx[seq(1, length(cx), by = 3)], fixed = T)
  Scx <- cx[seq(2, length(cx), by = 3)]
  Dcx <- cx[seq(3, length(cx), by = 3)]
  if (includes_MFE) {
    Mcx <- unlist(strsplit(Dcx, "[[:space:]]+"))
    Dcx <- Mcx[seq(1, length(Mcx), by = 2)]
    Mcx <- Mcx[seq(2, length(Mcx), by = 2)]
  }
  if (return_type == "DF") {
    return(data.frame("ID" = Ncx, "Seq" = Scx, "DBN" = Dcx, "MFE" = Mcx))
  }
  if (return_type == "BS") {
    OutBS <- BStringSet(x = Scx)
    names(OutBS) <- Ncx
    OutBS@elementMetadata <- DataFrame("DBN" = Dcx, "MFE" = Mcx) # OutBS["NC_001653.2 Hepatitis delta virus, complete genome"]@elementMetadata$DBN
    return(OutBS)
  }
}

PopForna <- function(DBN, pop = F, returl = T) { # To view on IMG/MER site.
  bu <- p0("http://nibiru.tbi.univie.ac.at/forna/forna.html?id=url/name&sequence=")
  bu_wd <- p0(bu, (unname(as.char(DBN))), "&structure=", DBN@elementMetadata$DBN)
  if (pop) {
    browseURL(bu_wd)
  }
  if (returl) {
    return(bu_wd)
  }
}

MMseqsClstsReader <- function(inpt = "resultsDB_clu.tsv", clsnm = "rdrp_id") {
  # Oldname: mmseqclstsvReader
  clsts <- fread(input = inpt, col.names = c("reprs", clsnm), data.table = F, skip = 0, header = F)
  singlts <- clsts[which(clsts$reprs == clsts[, clsnm]), ]
  clsts <- clsts[-which(clsts$reprs == clsts[, clsnm]), ]
  singlts <- singlts[-which(singlts$reprs %in% unique(clsts$reprs)), ]
  memtbl <- as.data.table(clsts)[, p0((.SD)), by = .(reprs)]
  memtbl$V1 <- gsub(x = memtbl$V1, pattern = 'c("', fixed = T, replacement = "")
  memtbl$V1 <- gsub(x = memtbl$V1, pattern = '")', fixed = T, replacement = "")
  memtbl$V1 <- gsub(x = memtbl$V1, pattern = '"', fixed = T, replacement = "")
  memtbl$V1 <- gsub(x = memtbl$V1, pattern = ")", fixed = T, replacement = "")
  memtbl$V1 <- gsub("\n", "", memtbl$V1, fixed = T)
  memtbl$V1 <- p0(memtbl$reprs, ", ", memtbl$V1)
  memtbl <- Rename1Col(memtbl, "V1", clsnm)
  singlts <- Rename1Col(singlts, "rdrp_id", clsnm)
  memtbl <- as.data.frame(rbind(memtbl, singlts))
  tstdf <- data.frame(do.call("rbind", strsplit(as.character(memtbl[, clsnm]), ", ", fixed = T)))
  memtbl$memlist <- apply(tstdf, 1, unique)
  memtbl[, clsnm] <- NULL
  memtbl$Nseq <- as.integer(lapply((memtbl$memlist), length))
  memtbl <- Rename1Col(memtbl, "memlist", clsnm)
  return(distinct(memtbl))
}

ScreenLastz <- function(reference, querry) { # Forked from Jerome Ambroise script: https://github.com/JeromeAmbroise/Pathogenomics2/blob/master/R/screenfunction.R

  rando <- floor(runif(1, min = 1, max = 10000000))
  try(unlink(p0("temp.", rando), recursive = TRUE))
  dir.create(p0("temp.", rando), showWarnings = F)
  myarg <- paste0(reference, "[multiple] ", querry, " --ambiguous=iupac --notransition --strand=both --step=100 --nogapped ‑‑format=rdotplot > ", p0("temp.", rando), "/result.maf")
  system2(command = "lastz", args = myarg)

  last <- try(read.table(p0("temp.", rando, "/", "result.maf")), silent = T)
  if (class(last) == "data.frame") {
    start.stop <- as.numeric(na.omit(suppressWarnings(as.numeric(as.character(last$V1)))))
    start <- start.stop[seq(1, length(start.stop), by = 2)]
    stop <- start.stop[seq(2, length(start.stop), by = 2)]
    GR <- GenomicRanges::GRanges(seqnames = "seq", ranges = IRanges(start = start, end = stop))
    GR.disjoin <- disjoin(GR)

    hitlength <- sum(width(GR.disjoin))
    seqlength <- sum(width(readDNAStringSet(reference)))

    percentage <- 100 * round(hitlength / seqlength, 3)
  } else {
    (percentage <- 0)
  }
  try(unlink(p0("temp.", rando), recursive = TRUE))
  return(percentage)
}

Kalign3 <- function(xs, param = NULL) {
  wd <- tempdir()
  dir <- getwd()
  temp_file <- basename(tempfile(tmpdir = wd))
  on.exit({
    file.remove(Sys.glob(paste(temp_file, ".*", sep = "")))
    setwd(dir)
  })
  setwd(wd)
  infile <- p0(temp_file, ".in")
  outfile <- p0(temp_file, ".aln")
  writeXStringSet(xs, infile, append = FALSE, format = "fasta")
  system(p0("kalign3 -in ", infile, " -out ", outfile, " -f fasta"))
  if (is(xs, "DNAStringSet")) {
    r <- readDNAMultipleAlignment(outfile, format = "fasta")
  }
  if (is(xs, "RNAStringSet")) {
    r <- readRNAMultipleAlignment(outfile, format = "fasta")
  }
  if (is(xs, "AAStringSet")) {
    r <- readAAMultipleAlignment(outfile, format = "fasta")
  }
  return(r)
}


GenericRunner <- function(command = "echo", param = list("threads" = 4, "mem" = 120), RunInTmp = F, FlgType = "-", AsType = " ", ...) {
  if (RunInTmp) {
    origdir <- getwd()
    wd <- tempdir()
  }
  if (!RunInTmp) {
    origdir <- getwd()
    wd <- origdir
  }
  setwd(wd)
  temp_file <- basename(tempfile(tmpdir = wd))
  on.exit({
    file.remove(Sys.glob(paste(temp_file, ".*", sep = "")))
    setwd(origdir)
  })
  Comsdf <- data.frame("argname" = names((param)), "argval" = unname(unlist(param)))
  Comsdf$astr <- paste0(FlgType, Comsdf$argname)
  Comsdf$astr <- paste0(Comsdf$astr, AsType, Comsdf$argval)
  Comscaf <- paste(sep = " ", command, paste(Comsdf$astr, collapse = " "))
  system(Comscaf, ...)
}

GenericRunnerInOut <- function(xs, infile = "infile.faa", inflag = "i", outflag = "o", outfile = "outfile.afa", keepinout = F, command = "echo", param = list("threads" = 4, "mem" = 120), RunInTmp = F, FlgType = "-", AsType = " ", ...) {
  # For biosequences.
  writeXStringSet(xs, infile, append = FALSE, format = "fasta")
  iol <- list(infile, outfile)
  names(iol) <- c(inflag, outflag)
  param <- append(param, iol)
  GenericRunner(param = param, ...)
  if (is(xs, "DNAStringSet")) {
    r <- readDNAMultipleAlignment(outfile, format = "fasta")
  }
  if (is(xs, "RNAStringSet")) {
    r <- readRNAMultipleAlignment(outfile, format = "fasta")
  }
  if (is(xs, "AAStringSet")) {
    r <- readAAMultipleAlignment(outfile, format = "fasta")
  }
  return(r)
}

###### ML/lm related  ######
GenerateSeqDescriptors <- function(indf = allmots, PrecentRand = 10, aacol = "AA_seq", HEC_Prob = F, AA_Prob = F, K.HEC = 3, K.aa = 2, casecol = "case", Exmp_src_faa = rdrps_faa, SEED = 123, valuecol = "motif_type", RemovesSparseCols = T) {
  # Oldname: GenerateMotifDescriptors
  # Need some refactoring but might work.
  set.seed(SEED)
  # indf=
  setDT(indf)
  trainR <- which(indf[, ..casecol] == "Train")
  testR <- setdiff(c(1:nrow(indf)), trainR)
  if (PrecentRand != 0) {
    Nfakes <- (PrecentRand / 100) * length(trainR)
    fake_seqs <- universalmotif::shuffle_sequences(sequences = Exmp_src_faa[(sample(x = length(Exmp_src_faa), size = Nfakes))])
    fakemotifs <- narrow(fake_seqs, start = 1, end = (sample(x = indf$mL[trainR], size = Nfakes)))
    names(fakemotifs) <- c(p0("fake.", (1:Nfakes)))
    fakemotifs <- XString2DF(fakemotifs, input_was = "new_name")
    fakemotifs[, valuecol] <- "666"
    fakemotifs$name_seq <- p0(fakemotifs$new_name, ".", unlist(fakemotifs[, valuecol]))
    fakemotifs <- Rename1Col(fakemotifs, "seq", aacol)
    fakemotifs <- Rename1Col(fakemotifs, "Length", "mL")
    fakemotifs$case <- "Train"
    indf <- rbind(indf, fakemotifs)
    trainR <- union(trainR, grep("fake", indf$new_name, fixed = T))
  }
  indf$pssp <- DECIPHER::PredictHEC(AAStringSet(unlist(indf[, ..aacol])))
  pssp <- BStringSet(indf$pssp)
  names(pssp) <- indf$name_seq
  psspdfdt <- XString2KmerDT(infaa = pssp, abcde = c("H", "E", "C"), k = K.HEC, out_probablity = HEC_Prob, input_was = "name_seq")
  colnames(psspdfdt) <- tolower(colnames(psspdfdt))
  AAset <- BStringSet(AAStringSet(unlist(indf[, ..aacol])))
  names(AAset) <- indf$name_seq
  AAdfdt <- XString2KmerDT(infaa = AAset, abcde = setdiff(AA_ALPHABET, c(".", "+", "-", "*")), k = K.aa, out_probablity = AA_Prob, input_was = "name_seq")
  indf <- distinct(merge(merge(indf, AAdfdt, by = "name_seq", all.x = T, all.y = F), psspdfdt, by = "name_seq", all.x = T, all.y = F))
  names(indf) <- trimws(names(indf)) # Legacy
  if (RemovesSparseCols == T) {
    indf <- setDT(RemovesSparseCols(as.data.frame(indf)))
  }
  vals <- unique(unlist(indf[, ..valuecol]))
  indf[, as.character(c(vals))] <- 0
  for (val in vals) {
    set(indf, i = which((indf[, ..valuecol] == val)), j = val, value = 1)
  }
  gc()
  return(indf)
}
