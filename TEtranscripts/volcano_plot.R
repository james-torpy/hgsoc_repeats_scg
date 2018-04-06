### 13.DE_FT_vs_HGSOC.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:

# Run on cluster with:
#briR
#qsub -N EDAEdgeR -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/DE \
#-j y -R y -pe smp 2 -V "Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/repeatsPipeline/exp9/13.DE_master.R"

### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "TEtranscripts"

descrip <- "TEtranscripts_primary_HGSOC_vs_FT"

# specify what combination of repeat genes (repeats), epigenetic modulators (epiMods),
# RNAi genes (RNAi) and protein-coding genes (pCoding) should contribute to the results:
#resultTypes <- c("repeats", "epiMods")
resultTypes <- c("repeats", "all")

# specify what padj and log2 fold change thresholds to use:
padjthresh <- NA
FCthresh <- 1

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")

# specify epigenetic modifier genes to include if necessary:
epiIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305", "ENSG00000276043", 
            "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", "ENSG00000101945")

epiSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", "SUV39H1")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0(projectDir, "/RNA-seq/raw_files/fullsamples/bowtell_primary/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/DEverify/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/DEverify/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


### 1. Load in differential expression results and separate repeats and other genes ###

allGenes <- read.table(paste0(resultsDir, "/TEtranscripts/test1/TEtranscripts_out_gene_TE_analysis.txt"),
                       header=T)
rownames(allGenes) <- allGenes$id
allGenes <- dplyr::select(allGenes, -id)

if ("repeats" %in% resultTypes) {
  # define repeat and sig DE repeat dfs:
  repGenes <- allGenes[grep("ENS",  rownames(allGenes), invert = T),]
  print(repGenes)
  
  if ( is.na(FCthresh) ) {
    sigGenes <- filter(repGenes, padj < padjthresh)
    repGenes$threshold <- as.factor(repGenes$padj < padjthresh)
  } else if ( is.na(padjthresh) ) {
    sigGenes <- repGenes[(repGenes$log2FoldChange > FCthresh)|(repGenes$log2FoldChange < -(FCthresh)), ]
    repGenes$threshold <- as.factor( (repGenes$log2FoldChange > FCthresh)|(repGenes$log2FoldChange < -(FCthresh)) )
  } else {
    sigGenes <- filter(repGenes, (padj < padjthresh & log2FoldChange < -(FCthresh))|(padj < padjthresh & log2FoldChange 
                                                                                     > FCthresh))
    repGenes$threshold <- as.factor((repGenes$padj < padjthresh & repGenes$log2FoldChange < -(FCthresh))|(repGenes$padj 
                                                                    <  padjthresh & repGenes$log2FoldChange > FCthresh))
  }
  
  ######
  #sigGenes <- filter(repGenes, (padj < padjthresh | log2FoldChange < -(FCthresh))|(padj < padjthresh | log2FoldChange >                                                                                   FCthresh))
  #repGenes$threshold <- as.factor((repGenes$padj < padjthresh | repGenes$log2FoldChange < -(FCthresh))|(repGenes$padj <  
  #                                                                   padjthresh | repGenes$log2FoldChange > FCthresh))
  ######
  
  sig <- subset(repGenes, threshold == T)
  
  # include the control genes for labelling:
  for (j in 1:length(posGeneIDs)) {
    if (j==1) {
      posGenes <- allGenes[ posGeneIDs[j],]
    } else {
      posGenes <- rbind(posGenes,   allGenes[posGeneIDs[j],])
    }
  }
  rownames(posGenes) <- posGeneNames
  
  for (j in 1:length(negGeneIDs)) {
    if (j==1) {
      negGenes <- allGenes[ negGeneIDs[j],]
    } else {
      negGenes <- rbind(negGenes,   allGenes[negGeneIDs[j],])
    }
  }
  rownames(negGenes) <- negGeneNames
  
  
  # set default threshold statuses  for control genes:
  ctlPadjthresh <- 0.05
  
  posGenes$threshold <- "POSITIVE"
  if (nrow(posGenes[posGenes$padj< ctlPadjthresh,])>0) {
    posGenes[posGenes$padj < ctlPadjthresh,]$threshold <- "POSSIG"
  }
  
  negGenes$threshold = "NEGATIVE"
  if (nrow(negGenes[negGenes$padj< ctlPadjthresh,])>0) {
    negGenes[negGenes$padj <ctlPadjthresh,]$threshold <-  "NEGSIG"
  }
  
  # load in interesting repeats, DE in primary HGSOC vs FT:
  intReps <- readRDS(file="/Users/jamestorpy/clusterHome//projects/hgsoc_repeats/RNA-seq/Robjects/exp9//sigReps.rds")
  intReps <- rownames(intReps)
  
  lab <- rbind(rbind(sig,   posGenes), negGenes)
  repGenes <- rbind(rbind(repGenes,   posGenes), negGenes)
  lab$genes <- rownames(lab)
  lab <- lab[is.finite(lab$log2FoldChange),]
  lab <- lab[lab$genes %in% intReps,]
  
  # plot on volcano plot:
  p <- ggplot(data=repGenes, aes( x=log2FoldChange, y=-log10(padj),    color=threshold))
  p <- p + geom_point(data=repGenes)
  p <- p + geom_text_repel(data=lab, aes(label=genes))
  p <- p + theme(legend.position =  "none")
  p <- p + labs(x="log2 fold change   vs FT control", y="-log10   padj")
  #p <- p +  xlim(c(-4, 4))
      
  if (length(FCthresh) == 0) {
    if (file.exists(paste0(plotDir,   "/", "_volcano_padj_",   padjthresh, "_", ".pdf"))) {
      print(paste0(plotDir, "/",  "_volcano_padj_",   padjthresh, "_", ".pdf"))
      p
    } else {
      print(paste0("Creating  ",plotDir, "/",   "_volcano_padj_", padjthresh, "_", ".pdf"))
      pdf(file = paste0(plotDir, "/",   "_volcano_padj_",  padjthresh, "_", ".pdf"))
      print(p)
      dev.off()
    }
  } else {
    if (file.exists(paste0(plotDir, "/",  "_volcano_padj",   padjthresh, "_FC", FCthresh, "_", ".pdf"))) {
      print(paste0(plotDir, "/",  "_volcano_padj",   padjthresh, "_FC", FCthresh, "_", ".pdf already exists"))
      p
    } else {
      print(paste0("Creating  ", plotDir, "/",  "_volcano_padj", padjthresh, "_FC", FCthresh, "_", ".pdf"))
      pdf(file = paste0(plotDir, "/",  "_volcano_padj",   padjthresh, "_FC", FCthresh, "_", ".pdf"))
      print(p)
      dev.off()
    }
  }
}    

if ("all" %in% resultTypes) {
  
  rownames(allGenes) <- gsub("\\..*$", "", rownames(allGenes))
  
  if (length(FCthresh) == 0) {
    sigGenes <- allGenes %>%
      rownames_to_column() %>%
      filter(padj < padjthresh) %>%
      column_to_rownames()
    allGenes$threshold <- as.factor(allGenes$padj < padjthresh)
  } else {
    sigGenes <- filter(allGenes, (padj < padjthresh & log2FoldChange < -(FCthresh))|(padj < padjthresh & log2FoldChange > FCthresh))
    allGenes$threshold <- as.factor((allGenes$padj < padjthresh & allGenes$log2FoldChange < -(FCthresh))|(allGenes$padj <  padjthresh & allGenes$log2FoldChange > FCthresh))
  }
  
  sig <- subset(allGenes, threshold == T)
  
  # include the control genes for labelling:
  for (j in 1:length(posGeneIDs)) {
    if (j==1) {
      posGenes <- allGenes[ posGeneIDs[j],]
    } else {
      posGenes <- rbind(posGenes,   allGenes[posGeneIDs[j],])
    }
  }
  rownames(posGenes) <- posGeneNames
  
  for (j in 1:length(negGeneIDs)) {
    if (j==1) {
      negGenes <- allGenes[ negGeneIDs[j],]
    } else {
      negGenes <- rbind(negGenes,   allGenes[negGeneIDs[j],])
    }
  }
  rownames(negGenes) <- negGeneNames
  
  # set default threshold statuses  for control genes:
  posGenes$threshold <- "POSITIVE"
  if (nrow(posGenes[posGenes$padj< padjthresh,])>0) {
    posGenes[posGenes$padj<  padjthresh,]$threshold <- "POSSIG"
  }
  
  negGenes$threshold = "NEGATIVE"
  if (nrow(negGenes[negGenes$padj< padjthresh,])>0) {
    negGenes[negGenes$padj<  padjthresh,]$threshold <-  "NEGSIG"
  }
  
  lab <- rbind(rbind(sig, posGenes), negGenes)
  lab <- lab[( lab$log2FoldChange > 10 | lab$log2FoldChange < -10 | lab$padj < 10e-15 ),]
  allGenes <- rbind(rbind(allGenes,   posGenes), negGenes)
  lab$genes <- rownames(lab)
  
  # add gene symbol annotations where relevant:
  egENSEMBL <- toTable(org.Hs.egENSEMBL)
  egSYMBOL <- toTable(org.Hs.egSYMBOL)
  
  # for rows with ensembl ids, annotate entrez ids and symbols in separate columns
  # for lab and sig:
  lab$gene_id <- egENSEMBL$gene_id[match(rownames(lab), egENSEMBL$ensembl_id)]
  lab$symbol <- egSYMBOL$symbol[match(lab$gene_id, egSYMBOL$gene_id)]
  
  sig$gene_id <- egENSEMBL$gene_id[match(rownames(sig), egENSEMBL$ensembl_id)]
  sig$symbol <- egSYMBOL$symbol[match(sig$gene_id, egSYMBOL$gene_id)]
  
  allGenes$gene_id <- egENSEMBL$gene_id[match(rownames(allGenes), egENSEMBL$ensembl_id)]
  allGenes$symbol <- egSYMBOL$symbol[match(allGenes$gene_id, egSYMBOL$gene_id)]
  
  # plot on volcano plot:
  p <- ggplot(data=allGenes, aes( x=log2FoldChange, y=-log10(padj), color=threshold) )
  p <- p + geom_point(data=allGenes)
  p <- p + geom_text_repel(data=lab, aes(label=symbol))
  p <- p + theme(legend.position =  "none")
  p <- p + labs(x="log2 fold change   vs FT control", y="-log10   padj")
  if (length(FCthresh) == 0) {
    if (file.exists(paste0(plotDir,   "/", "_volcano_padj_10e_neg15_", "_allGenes.pdf"))) {
      print(paste0(plotDir, "/",  "_volcano_padj_10e_neg15_", "_allGenes.pdf"))
      p
    } else {
      print(paste0("Creating  ",plotDir, "/",   "_volcano_padj_10e_neg15_", "_allGenes.pdf"))
      pdf(file = paste0(plotDir, "/",   "_volcano_padj_10e_neg15_", "_allGenes.pdf"))
      print(p)
      dev.off()
    }
  } else {
    if (file.exists(paste0(plotDir, "/",  "_volcano_padj_10e_neg15_", "_allGenes.pdf"))) {
      print(paste0(plotDir, "/",  "_volcano_padj_10e_neg15_", "_allGenes.pdf already exists"))
      p
    } else {
      print(paste0("Creating  ", plotDir, "/",  "_volcano_padj_10e_neg15_", "_allGenes.pdf"))
      pdf(file = paste0(plotDir, "/",  "_volcano_padj_10e_neg15_", "_allGenes.pdf"))
      print(p)
      dev.off()
    }
  }
}
