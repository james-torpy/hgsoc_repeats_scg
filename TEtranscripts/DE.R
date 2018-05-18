### DE.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:


### 0. Define variables/paths ###

# load packages needed:
library(tibble)
library(dplyr)
library(RColorBrewer)
library(RUVSeq)
library(ggplot2)
library(ggrepel)
library(preprocessCore)
library(edgeR)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "TEtranscripts"
Type <- "custom3"

diff_exp <- "primary_HGSOC_vs_FT"
descrip <- "TEtranscripts_EdgeR_primary_HGSOC_vs_FT"
count_tool <- "TETranscripts"

# specify what combination of repeat genes (repeats) and other genes,
# (all, both, other) should contribute to the results:
resultTypes <- c("repeats", "all", "both", "other")

# specify what FDR and log2 fold change thresholds to use:
Pthresh <- 0.001
FCthresh <- 2

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")

# specify other genes to include if necessary:
#otherIDs <- c("ENSG00000130816", "ENSG00000119772", "ENSG00000088305", 
#	"ENSG00000276043", "ENSG00000138336", "ENSG00000168769", "ENSG00000187605", 
# "ENSG00000101945")

#otherSym <- c("DNMT1", "DNMT3A", "DNMT3B", "UHRF1", "TET1", "TET2", "TET3", 
#	"SUV39H1")

#otherIDs <- c("ENSG00000128383", "ENSG00000179750", "ENSG00000128394")
#otherSym <- c("APOBEC3A", "APOBEC3B", "APOBEC3F")

#otherIDs <- c("ENSG00000104824", " ENSG00000136436", "ENSG00000161011", 
#	"ENSG00000125207", "ENSG00000197181", "ENSG00000184571", "ENSG00000134627")
#otherSym <- c("HNRNPL", "NPD52", "p62", "PIWIL1", 
#              "PIWIL2", "PIWIL3", "PIWIL4")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
rawDir <- paste0(projectDir, 
                 "/RNA-seq/raw_files/fullsamples/bowtell_primary/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
inDir <- paste0(resultsDir, "/", expName, "/", diff_exp, "/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


################################################################################
### 1. Load in DE results and sort into all genes and repeats ###
################################################################################

allGenes <- read.table(paste0(inDir, "/primary_HGSOC_vs_FT_sigdiff_gene_TE.txt"), header = T)
allGenes <- subset(allGenes, select = c(id, log2FoldChange, padj))
colnames(allGenes) <- c("gene_id", "logFC", "padj")
rownames(allGenes) <- gsub("\\..*$", "", allGenes$gene_id)

        
# annotate allGenes with entrez ids and symbols in separate columns:
egENSEMBL <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

allGenes$gene_id <- egENSEMBL$gene_id[match(rownames(allGenes), egENSEMBL$ensembl_id)]
allGenes$symbol <- egSYMBOL$symbol[match(allGenes$gene_id, egSYMBOL$gene_id)]

# define repeat and sig DE repeat dfs:
repGenes <- allGenes[grep("ENS",  rownames(allGenes), invert = T),]
repGenes$symbol <- rownames(repGenes)
repGenes <- subset(repGenes, select = -gene_id)
repGenes$threshold <- NA
repGenes$threshold[grep("LINE|SINE|ALU", repGenes$symbol)] <- "label"
repGenes$threshold[repGenes$logFC > 0] <- "label"
repGenes$threshold[repGenes$padj < Pthresh & repGenes$logFC < -(FCthresh)] <- "label"
print(repGenes)

labGenes <- repGenes[!is.na(repGenes$threshold),]

################################################################################
### 2. Create volcano plots ###
################################################################################

p <- ggplot(data=repGenes, aes(x=logFC, y=-log10(padj)))
p <- p + geom_point(data=repGenes)
p <- p + geom_text_repel(data=labGenes, aes(label=symbol))
p <- p + theme(legend.position = "none")
p <- p + labs(x="log2 fold change vs FT control", y="-log10 adjusted p-value")
#p <- p +  xlim(c(-5, 5))
p





