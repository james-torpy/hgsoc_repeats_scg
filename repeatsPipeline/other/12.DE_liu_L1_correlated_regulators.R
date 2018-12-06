### 13.DE_master.R ###

# This script takes a list of dfs of different classes and types of repeat counts for
# HGSOC and FT control RNA-seq data sets and performs DE analysis:

# Run on cluster with:
#briR
#qsub -N salDE -b y -wd \
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
library(preprocessCore)
library(edgeR)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
Type <- "custom3"


################################################################################
### Options ###
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC_vs_FT ###
# sTypes <- c("FT", "HGSOC")
# sGroups <- list("FT", c("prPT", "rfPT", "arPT", "mrPT", "erPT"))
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_primary_HGSOC_vs_FT_with_correlated_liu_L1_regulators"
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC_vs_HRD ###

# sTypes <- c("HRD", "CCNEamp", "unknown")
# descrip <- "htseq_EdgeR_primary_HGSOC_vs_HRD"
# sGroups <- list("HRD", "CCNEamp", "unknown")
# names(sGroups) <- sTypes
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC__unknown_driver_vs_known_driver_with_liu_L1_regulators ###

# sTypes <- c("known", "unknown")
# descrip <- "htseq_EdgeR_primary_HGSOC_unknown_driver_vs_known_driver_with_liu_L1_regulators"
# sGroups <- list(c("bothDrivers", "CCNEamp", "HRD"), "unknown")
# names(sGroups) <- sTypes
################################################################################

################################################################################
### htseq_EdgeR_HGSOC_drug_cats_vs_FT ###

# sTypes <- c("FT", "primary_resistant", "acquired_resistant", "drug_responders",
#  "recurrent_ascites", "metastatic")
# sGroups <- list("FT", "prPT", "rfPT", "arPT", "mrPT", "erPT", "rcAF", "msST")
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_HGSOC_drug_cats_vs_FT"
################################################################################

################################################################################
### htseq_EdgeR_HGSOC_prim_drug_cats ###

sTypes <- c("primary_resistant", "acquired_resistant", "multiple_responders",
            "extreme_responders")
sGroups <- list(c("prPT", "rfPT"), "arPT", "mrPT", "erPT")
names(sGroups) <- sTypes
descrip <- "htseq_EdgeR_HGSOC_prim_drug_cats"
################################################################################

################################################################################
### htseq_EdgeR_HGSOC_resistant_vs_responsive ###

# sTypes <- c("acquired_resistant", "primary_resistant", "drug_responders")
# sGroups <- list("arPT", c("prPT", "rfPT"), c("mrPT", "erPT"))
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_HGSOC_resistant_vs_responsive"
################################################################################

################################################################################
### SalmonTE_primary_HGSOC_CCNEamp_vs_HRD ###

#sTypes <- c("bothDrivers", "FT", "HRD", "CCNEamp", "unknown_driver")
#descrip <- "SalmonTE_primary_HGSOC_CCNEamp_vs_HRD"
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC_vs_FT_with_liu_L1_regulators ###
# sTypes <- c("FT", "HGSOC")
# sGroups <- list("FT", c("prPT", "rfPT", "arPT", "mrPT", "erPT"))
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_primary_HGSOC_vs_FT_with_liu_L1_regulators"
# otherType <- c("suppressor", "activator")
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC_vs_CCNE1_with_liu_L1_regulators ###
# sTypes <- c("CCNEamp", "HRD", "unknown")
# sGroups <- list("CCNEamp", "HRD", "unknown")
# names(sGroups) <- sTypes
# descrip <- "htseq_EdgeR_primary_HGSOC_vs_CCNE1_with_liu_L1_regulators"
# otherType <- c("suppressor", "activator")
################################################################################

################################################################################
### htseq_EdgeR_primary-resistant_HGSOC_vs_non-primary-resistant_HGSOC_with_liu_L1_regulators ###
# sTypes <- c("non-pr_HGSOC", "pr_HGSOC")
# sGroups <- list(c("arPT", "mrPT", "erPT"), c("prPT", "rfPT"))
# names(sGroups) <- sTypes
# descrip <-
#   "htseq_EdgeR_primary-resistant_HGSOC_vs_non-primary-resistant_HGSOC_with_liu_L1_regulators"
################################################################################

# define comparison parameters:
primaryOnly <- TRUE
cat_by_driver <- FALSE
EDAnormalise <- FALSE
count_tool <- "EdgeR"
customSamples <- FALSE
include_ctls <- FALSE

# define custom samples if needed:
#cus <- c("")

# specify what combination of repeat genes (repeats) and other genes,
# (all, both, other) should contribute to the results:
resultTypes <- c("repeats", "other")

# define sample group to use as control:
ctl <- "acquired_resistant"

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.3
FCthresh <- 0

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
rawDir <- paste0(projectDir, 
                 "/RNA-seq/raw_files/fullsamples/bowtell_primary/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))

# specify other genes to include if necessary:
other_df <- read.table(paste0(refDir, "/correlated_liu_2018_L1_regulators.txt"), 
  header=T, sep = "\t")

otherIDs <- other_df$ensembl_id
otherSym <- other_df$symbol
suppressor_df <- other_df[other_df$effect == "suppressor",]
activator_df <- other_df[other_df$effect == "activator",]


################################################################################
### 1. Load in all counts ###
################################################################################

# define functions for this section:
counts_bind <- function(counts1, counts2) {
  # append counts1 to counts2:
  counts_all <- rbind(custom3Counts, gcCounts)
  
  # make rownames gene_id, get rid of latter column and change
  # storage mode from factor to integer:
  rownames(counts_all) <- counts_all$gene_id
  return(subset(counts_all, select=-gene_id))
}

if ( count_tool=="EdgeR" ) {
  
  if ( !file.exists(paste0(RobjectDir, "/", Type, "_", count_tool, 
                           "_counts.rds")) ) {
    
    writeLines("\n")
    print("EdgeR counts data frame does not exist, creating now...")
    
    custom3Counts <- readRDS(paste0(RobjectDir, "/", Type, 
                                    "_allcounts.htseq.rds"))
    gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
    
    # append gcCounts to custom3Counts:
    Counts <- counts_bind(custom3Counts, gcCounts)
    
    saveRDS(Counts, paste0(RobjectDir, "/", Type, "_", count_tool, 
                           "_counts.rds"))
  } else {
    
    print("Loading EdgeR counts data frame...")
    Counts <- readRDS(paste0(RobjectDir, "/", Type, "_", count_tool, 
                           "_counts.rds"))
  }
  
} else if ( count_tool=="SalmonTE" ) {
  
  if ( !file.exists(paste0(RobjectDir, "/", Type, 
                           "_", count_tool, "_counts.rds")) ) {
    
    writeLines("\n")
    print("SalmonTE counts data frame does not exist, creating now...")
    
    custom3Counts <- readRDS(paste0(RobjectDir, "/", Type, 
                                    "_counts.SalmonTE.rds"))
    gcCounts <- readRDS(paste0(RobjectDir, "/gc_counts.Salmon.rds"))
    
    # append gcCounts to custom3Counts:
    Counts <- counts_bind(custom3Counts, gcCounts)
    
    saveRDS(Counts, (paste0(RobjectDir, "/", Type, "SalmonTE_counts.rds")))
    
  } else {
    
    print("Loading SalmonTE counts data frame...")
    Counts <- readRDS(file=paste0(RobjectDir, "/", Type, 
                                  "SalmonTE_counts.rds"))
  }
}
    
# if necessary, select primary samples only:
if ( primaryOnly == TRUE) {
    Counts <- Counts[,grep("PT|FT", colnames(Counts))]
  }

# select custom samples:
if (customSamples) {
  Counts <- Counts[,colnames(Counts) %in% cus]
}

# re-categorize samples as HRD, CCNE_amp, both_drivers or unknown_driver:
if ( cat_by_driver ) {
  # load in sample key for categories homologous repair deficient (HRD) and 
  #cyclin E gain/amplification (CCNE):
  HRDkey <- read.table(paste0(rawDir, "/hrd_samples.txt"), header=F, sep="\t")
  HRDnos <- gsub("AOCS_", "", HRDkey[,1])
  CCNEkey <- read.table(paste0(rawDir, "/ccne_gain_or_amp_samples.txt"), 
                        header=F, sep="\t")
  CCNEnos <- gsub("AOCS_", "", CCNEkey[,1])
  Counts <- Counts[,grep("rcAF|pAF|msST", colnames(Counts), invert=T)]
  for (i in 1:length(colnames(Counts))) {
    print(i)
    if (length(grep("FT", colnames(Counts)[i]))<1) {
      no <- gsub(
        "_[a-zA-Z].*$", "", gsub("AOCS_", "", colnames(Counts)[i])
      )
      print(paste0("ID number is ", no))
      if (no %in% HRDnos & no %in% CCNEnos) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_bothDrivers")
      } else if (no %in% HRDnos & !(no %in% CCNEnos)) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_HRD")
      } else if (no %in% CCNEnos & !(no %in% HRDnos)) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_CCNEamp")
      } else if (!(no %in% HRDnos | no %in% CCNEnos)) {
        colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_unknown")
      }
      colnames(Counts)[i] <- gsub("[a-z].*[A-Z][A-Z]_", "", 
                                  colnames(Counts)[i])
    }
  }
} else {
  
  # remove any samples not belonging to any group:
  Counts <- Counts[, colnames(
    Counts[, gsub(
      ".*\\_", "", colnames(Counts)
    ) %in% unlist(sGroups)]
  )]
}

# change sample names according to grouping:
for (i in 1:length(sGroups)) {
  for (n in sGroups[[i]]) {
    colnames(Counts) <- gsub(n, names(sGroups)[i], colnames(Counts))
  }
}

# remove samples not belonging to any group:
if ( cat_by_driver ) {
  Counts <- Counts[,gsub("AOCS_[0-9][0-9][0-9]_", "", colnames(Counts)) %in% names(sGroups)]
} else {
  Counts <- Counts[,gsub("AOCS_[0-9][0-9][0-9]_", "", colnames(Counts)) %in% names(sGroups)]
}
  
# eliminate lowly expressed genes (rows where there are less than 3 counts 
# where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= (ncol(Counts)/3)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))
  
  
############################################################################
### 2. Perform pre-normalisation PCA and RLE plots ###
############################################################################
  
# create pre-normalised PCA plot from counts and plot:
if (ncol(Counts) > nrow(Counts)) {
  pca <- prcomp(Counts)
} else {
  pca <- princomp(Counts)	  
}
  
if (file.exists(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf already exists,
               no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaCompsPrenormGC.pdf"))
  plot(pca)
  dev.off()
}
  
# change the order of columns of Counts to alphabetical order:
Counts <- Counts[,order(
  gsub(
    "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
  )
)]
  
# define sample groups:
splt <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
      )
    ), length
  )
)
  
for (i in 1:length(splt)) {
  if (i==1) {
    typeF <- c(rep(names(splt)[i], splt[i]))
  } else {
    typeF <- c(typeF, rep(names(splt)[i], splt[i]))
  }
}
levels(typeF) <- sTypes
  
sampleNos <- unlist(
  lapply(
    split(
      colnames(Counts), gsub(
        "\\.1", "", gsub(
          "AOCS.*_[0-9][0-9][0-9]_", "", colnames(Counts)
        )
      )
    ), length
  )
)

saveRDS(sampleNos, file = paste0(newRobjectDir, "/sample_no_per_cat.rds"))

# convert Counts into SeqExpressionSet - elements need to be delisted and 
# changed to integers first:
Counts <- apply(Counts, 2, unlist)
storage.mode(Counts) <- "integer"
set <- newSeqExpressionSet(Counts, phenoData = data.frame(typeF, 
                                                          row.names=colnames(Counts)))

# create pre-norm RLE plot:
if (file.exists(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf already exists, 
             no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  par(mar=c(1,1,1,1))
  pdf(file = paste0(plotDir, "/", Type, "_RLEPrenormGC.pdf"))
  plotRLE(set)
  dev.off()
}

# create RUVseq pre-norm PCA:
if (file.exists(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"))) {
  print(paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf already exists, 
             no need to create"))
} else {
  print(paste0("Creating ", plotDir, "/", Type, "_pcaPrenormGC.pdf"))
  pdf(file = paste0(plotDir, "/", Type, "_pcaPrenormGC.pdf"), height = 10, 
      width = 12)
  plotPCA(set, cex=0.7)
  dev.off()
}
  

##############################################################################
### 3. perform normalisation and DE on counts:
##############################################################################

  
# perform between lane full normalisation:
y <- DGEList(counts = Counts, group = typeF)

# normalise for library size:
y <- calcNormFactors(y)


# create an MDS plot to show relative similarities of the samples and save to Dir:
if (file.exists(paste0(plotDir, "/edger_MDS.pdf"))) {
  paste0(plotDir, "/edger_MDS.pdf already exists")
  pdf(paste0(plotDir, "/edger_MDS.pdf"),width=16,height=12)
  plotMDS(y)
} else {
  paste0("Generating ", plotDir, "/edger_MDS.pdf")
  pdf(paste0(plotDir, "/edger_MDS.pdf"),width=16,height=12)
  plotMDS(y)
  dev.off()
}

# design matrix labelling all sample types:
design <- model.matrix(~0+typeF)
colnames(design) <- gsub("typeF", "", colnames(design))

# estimate dispersion:
disp <- estimateDisp(y, design=design)

# adjust values using dispersion:
fit <- glmFit(disp, design=design, robust=TRUE)

saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "DEfit.rds"))
save.image(paste0(newRobjectDir, "/", Type, "DEdone.rds"))

# determine which column has control:
ctlInd <- which(colnames(design)==ctl)
  
con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))

# put sTypes in alphabetical order:
sTypes <- sTypes[order(sTypes)]


# check parameters and give the user the option to continue or not:
writeLines("\n")
print("Contrast is: ")
print(con)
print("Column names of design are: ")
print(colnames(design))
print("Design matrix is: ")
print(design)

for (i in 1:ncol(design)) {
  print(i)
  if ( i!=ctlInd ) {
    comp <- paste0(sTypes[i], "_vs_", ctl)
    
    # perform likelihood ratio test:
    con[i] <- 1
    lrt <- glmLRT(fit, contrast = con)
    
    # determine the top DE genes:
    topTags(lrt)
    
    if (file.exists(paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))) {
      print(paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds already exists, 
                   no need to create"))
    } else {
      print(paste0("Creating ", newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))
      saveRDS(lrt, file = paste0(newRobjectDir, sTypes[i], "_vs_", ctl, "_lrt.rds"))
    }
    
    # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
    DEs <- summary(result <- decideTestsDGE((lrt)))
    
    # fetch all gene DE info, 
    allGenes <- as.data.frame(topTags(lrt, n=Inf))
    
    # annotate allGenes with entrez ids and symbols in separate columns:
    egENSEMBL <- toTable(org.Hs.egENSEMBL)
    egSYMBOL <- toTable(org.Hs.egSYMBOL)
    
    allGenes$gene_id <- egENSEMBL$gene_id[match(rownames(allGenes), egENSEMBL$ensembl_id)]
    allGenes$symbol <- egSYMBOL$symbol[match(allGenes$gene_id, egSYMBOL$gene_id)]

    if (!(ctlInd==1)) {
      if (i==1) {
        allGenesList <- list(allGenes)
      } else {
        allGenesList[[i]] <- allGenes
      }
      
    } else {
      if (i==2) {
        allGenesList <- list(allGenes)
      } else {
        allGenesList[[i]] <- allGenes
      }
    }

    # create threshold column for FC/FDR cutoff:
    if (length(FCthresh) == 0) {
      sigGenes <- filter(allGenes, FDR < FDRthresh)
      allGenes$threshold <- as.factor(allGenes$FDR < FDRthresh)
    } else {
      sigGenes <- filter(allGenes, (FDR < FDRthresh & logFC < -(FCthresh))|(FDR < FDRthresh & logFC > FCthresh))
      allGenes$threshold <- as.factor((allGenes$FDR < FDRthresh & allGenes$logFC < -(FCthresh))|(allGenes$FDR <  FDRthresh & allGenes$logFC > FCthresh))
    }
    
    
    ##############################################################################
    ### 4. Create DE data frames for repeats:
    ##############################################################################
    
    # define repeat and sig DE repeat dfs:
    repGenes <- allGenes[grep("ENS",  rownames(allGenes), invert = T),]
    print(repGenes)

    # add 'type' identifier column:
    repGenes$type <- "repeat"
    
    repGenes <- repGenes[grep("^L1", rownames(repGenes)),]
    repGenes$symbol <- rownames(repGenes)
    
    sig_rep <- subset(repGenes, threshold == T)
    
    if (!(ctlInd==1)) {
      if (i==1) {
        allReps <- list(repGenes)
      } else {
        allReps[[i]] <- repGenes
      }
      
      if (i==1) {
        sigReps <- list(sig_rep)
      } else {
        sigReps[[i]] <- sig_rep
      }
    } else {
      if (i==2) {
        allReps <- list(repGenes)
      } else {
        allReps[[i]] <- repGenes
      }
      
      if (i==2) {
        sigReps <- list(sig_rep)
      } else {
        sigReps[[i]] <- sig_rep
      }
    }
    
    
    ##############################################################################
    ### 5. Create DE data frames for gencode genes:
    ##############################################################################
    
    gcGenes <- allGenes[grep("ENS",  rownames(allGenes)),]

    # add 'type' identifier column:
    gcGenes$type <- "gc"

    sig_gc <- subset(allGenes, threshold == T)
    
    if (!(ctlInd==1)) {
      if (i==1) {
        sig_gc_GenesList <- list(sig_gc)
      } else {
        sig_gc_GenesList[[i]] <- sig_gc
      }
    } else {
      if (i==2) {
        sig_gc_GenesList <- list(sig_gc)
      } else {
        sig_gc_GenesList[[i]] <- sig_gc
      }
    }


    ##############################################################################
    ### 6. Create positive and negative control genes data frame:
    ##############################################################################

    # include the control genes for labelling:
    for (j in 1:length(posGeneIDs)) {
      if (j==1) {
        posGenes <- allGenes[ posGeneIDs[j],]
      } else {
        posGenes <- rbind(posGenes, allGenes[posGeneIDs[j],])
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
    
    # set default threshold statuses for control genes:
    posGenes$threshold <- "POSITIVE"
    if (nrow(posGenes[posGenes$FDR< FDRthresh,])>0) {
      posGenes[posGenes$FDR<  FDRthresh,]$threshold <- "POSSIG"
    }
    
    negGenes$threshold = "NEGATIVE"
    if (nrow(negGenes[negGenes$FDR< FDRthresh,])>0) {
      negGenes[negGenes$FDR<  FDRthresh,]$threshold <-  "NEGSIG"
    }

    ctlGenes <- rbind(posGenes, negGenes)

    # add 'type' identifier column:
    ctlGenes$type <- "ctl"


    ##############################################################################
    ### 7. Create volcano plots:
    ##############################################################################
        
    
    suppressor_df$ensembl_id <- as.character(suppressor_df$ensembl_id)
    suppressor_genes <- allGenes[suppressor_df$ensembl_id,]
    suppressor_genes$ensembl_id <- rownames(suppressor_genes)
    suppressor_genes <- merge(suppressor_df, suppressor_genes, by.x="ensembl_id", by.y="ensembl_id")
    suppressor_genes$type <- "suppressor"
    suppressor_genes <- subset(suppressor_genes, select=-symbol.y)
    
    activator_df$ensembl_id <- as.character(activator_df$ensembl_id)
    activator_genes <- allGenes[activator_df$ensembl_id,]
    activator_genes$ensembl_id <- rownames(activator_genes)
    activator_genes <- merge(activator_df, activator_genes, by.x="ensembl_id", by.y="ensembl_id")
    activator_genes$type <- "activator"
    activator_genes <- subset(activator_genes, select=-symbol.y)
    
    other_genes <- rbind(suppressor_genes, activator_genes)
    # combine 'threshold' and 'type' columns:
    other_genes$type_thresh <- paste0(other_genes$type, "_", other_genes$threshold)
    other_genes$type_thresh <- factor(other_genes$type_thresh, 
                                      levels = c("suppressor_FALSE", "activator_FALSE", "suppressor_TRUE", 
                                                 "activator_TRUE"))
    
    lab <- other_genes[other_genes$FDR < FDRthresh,]
    rownames(lab) <- lab$ensembl_id
    lab <- subset(lab, select = -c(ensembl_id, effect, type_thresh))
    lab <- lab[,c(2:7, 1, 8:9)]
    colnames(lab)[7] <- "symbol"
    lab <- rbind(lab, sig_rep)
        
    # combine 'threshold' and 'type' columns:
    repGenes$type_thresh <- paste0(repGenes$type, "_", repGenes$threshold)
    repGenes <- repGenes[grep("^L1", rownames(repGenes)),]
    other_genes <- subset(other_genes, select = -c(ensembl_id, effect))
    other_genes <- other_genes[,c(2:7, 1, 8:10)]
    colnames(other_genes)[7] <- "symbol"
    Genes <- rbind(repGenes, other_genes)
    
    lab$type_thresh <- paste0(lab$type, "_", lab$threshold)
    
    Genes$type_thresh <- factor(Genes$type_thresh, levels = c("repeat_FALSE", 
      "suppressor_FALSE", "activator_FALSE", "repeat_TRUE", "suppressor_TRUE",
      "activator_TRUE"))
    
    cols <- brewer.pal(6,"Set3")
    cols <- c(cols[1], cols)
    cols2 <- c("#B4D8B5", "lightpink", "#4CD64C", "deepskyblue", "#FB8072")
    
    # plot on volcano plot:
    p <- ggplot(data=Genes, aes( x=logFC, y=-log10(FDR), color=type_thresh))
    p <- p + geom_point(data=Genes)
    p <- p + geom_text_repel(data=lab, aes(label=symbol))
    # p <- p + theme(legend.position =  "none")
    p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
    p <- p + scale_colour_manual(values = cols2)
    p <- p +  xlim(c(-4, 4))
    
    pdf(file = paste0(plotDir, "/", comp, "_FDR_", FDRthresh, ".pdf"))  
    print(p)
    dev.off()
  }
  con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
}
