### 20.expression_vs_str_var_cats.R ###

# This script takes DE information for each primary HGSOC sampl vs pooled FT controls and plots it against number of SVs
# of different categories

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
descrip <- "expression_vs_structural_rearrangement_cats_indivHGSOC"
Type <- "custom3"
metric <- "HGSOC_vs_FT"
#repeatList <- c("(CATTC)n", "HSATII", "L1M3a", "L1M3c", "L1M4c", "L1M5", "L1M8", "L1MA5", "L1MC")
repeatPattern <- ""
FDRthresh <- 0.05
removeSamples <- c("pAF", "rcAF")

# define plot colours:
no <- 60
qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
rawDir <- paste0(projectDir, "/RNA-seq/raw_files")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
       expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))

# collect DE repeats between primary HGSOC and FT controls (EdgeR):
DE <- readRDS(file = paste0(RobjectDir, "/DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/custom3_DEsigReps.rds"))
repeatList <- rownames(DE[[1]])


### 1. Load in DESeq results and create logFC data frame ###

DElist <- readRDS(file=paste0(RobjectDir, "/htseq_single_HGSOCs_vs_FT_DElist.rds"))

# Concatenate all logFC values into data frame
i=1
logFC_list <- lapply(DElist, function(x) {
  df <- data.frame(x$log2FoldChange)
  rownames(df) <- x$id
  colnames(df) <- names(DElist)[i]
  df[df[,1]=="-Inf",] <- 0
  i <<- i+1
  return(df)
})

logFC <- do.call("cbind", logFC_list)

# only keep repeat information:
logFC <- logFC[grep("ENSG", rownames(logFC), invert = T),]

saveRDS(logFC, file=paste0(RobjectDir, "/repeat_logFC.rds"))

# Concatenate all padj values into data frame
i=1
padj_list <- lapply(DElist, function(x) {
  df <- data.frame(x$padj)
  rownames(df) <- x$id
  colnames(df) <- names(DElist)[i]
  i <<- i+1
  return(df)
})

padj <- do.call("cbind", padj_list)

# only keep repeat information:
padj <- logFC[grep("ENSG", rownames(padj), invert = T),]

saveRDS(logFC, file=paste0(RobjectDir, "/repeat_padj.rds"))

# Concatenate all FC values into data frame
i=1
FC_list <- lapply(DElist, function(x) {
  df <- data.frame(x$foldChange)
  rownames(df) <- x$id
  colnames(df) <- names(DElist)[i]
  df[df[,1]=="-Inf",] <- 0
  i <<- i+1
  return(df)
})

FC <- do.call("cbind", FC_list)

# only keep repeat information:
FC <- FC[grep("ENSG", rownames(FC), invert = T),]


### 2. Order samples by structural rearrangement number ###

# load structural variation counts:
str_var <- read.table(file=paste0(rawDir, "/fullsamples/bowtell_primary/SV_counts.txt"),
                      header = T, fill = T, sep = "\t")

# create a list with samplenames ordered by structural variation feature:
for (l in 1:length(colnames(str_var))) {
  print(l)
  if (l==1) {
    orders <- list(rownames(str_var[order(str_var[,l]),]))
    names(orders)[l] <- colnames(str_var)[l]
  } else {
    orders[[l]] <- rownames(str_var[order(str_var[,l]),])
    names(orders)[l] <- colnames(str_var)[l]
  }
}

for (m in 1:length(orders)) {
  
  # for each structural feature, order the logFC, FC and padj data frames by samples with
  # least to most of that feature:
  IDs <- paste0(orders[[m]], "_vs_FT")
  mtch <- na.omit(match(IDs, colnames(logFC)))
  logFC <- logFC[,mtch]
  padj <- padj[,mtch]
  FC <- FC[,mtch]
  
  for (j in 1:nrow(logFC)) {
    if (length(grep("ENSG", rownames(logFC)[j])) < 1) {
      rp <- as.data.frame(t(logFC[j,]))
      rp$sample <- factor(rownames(rp), levels = rownames(rp))
      rp$padj <- as.numeric(t(padj[j,]))
      rp$pThresh <- as.logical((rp$padj < 0.05))
      rp$sig <- NA
      rp$sig[rp$pThresh] <- "< 0.05"
      rp$sig[!(rp$pThresh)] <- "non-significant"
      
      rp$FC <- as.numeric(t(FC[j,]))
      rp$str_no <- str_var[gsub(
          "_vs_FT", "", rownames(rp),
        ), ][m]
      
      colnames(rp) <- c("log2FC", "sample", "padj", "pThresh", "sig", "FC", "str_no")
      if (j==1) {
        rpList <- list(rp)
      } else {
        rpList[[j]] <- rp
      }
    }
  }
  #rpList <- rpList[grep("ENSG", names(rpList), invert=T)]
  # reduce list to those DE in EdgeR primary HGSOC vs FT analysis:
  names(rpList) <- rownames(logFC)[1:length(rpList)]
  #rpList <- rpList[names(rpList) %in% repeatList]
  
  #saveRDS(rpList, file=paste0(RobjectDir, "/rpList.rds"))
  
  # check if repeat FC distribution is normal:
  #allFC <- lapply(rpList, function(x) {
  #  return(x$FC)
  #})
  #allFC <- unlist(allFC)
  #repFC <- allFC[grep("ENSG", names(allFC), invert=T)]
  
  #dFC <- density(repFC)
  #plot(dFC)
  
  n=1
  corList <- lapply(rpList, function(x) {
    df <- data.frame(x$str_no, x$FC)
    #df <- df[rpList[[n]]$sig == "< 0.05",]
    n <<- n+1
    return(cor(df, method = "spearman")[1,2])
  })
  names(corList) <- names(rpList)
  cors <- na.omit(unlist(corList))
  
  rpInt <- rpList[cors > 0.4 | cors < -0.4]
  
  
  ### 4. Plot logFC vs structural rearrangement for each repeat ###
  
  for (k in 1:length(rpInt)) {
    rpName <- gsub("\\(|\\)", "", names(rpInt)[k])
    pdf(file = paste0(plotDir, "/", rpName, "_", names(orders)[m], ".pdf"), width=10, height=10)
    p <- ggplot(rpInt[[k]], aes(x = factor(rpInt[[k]]$sample), y = rpInt[[k]]$log2FC))
    p <- p + geom_bar(stat="identity", position = "dodge", aes(fill = rpInt[[k]]$sig))
    p <- p + ylab("log2FC vs FT")
    p <- p + xlab("sample")
    p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
    p <- p + annotate("text", x=45, y=7, label=paste0("Spearman = ", as.character(cors[rpName])))
    p <- p + ggtitle(paste0(rpName, " vs ", names(orders)[m]))
    print(p)
    dev.off()
  }
  
}











