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
sTypes <- c("FT", "HGSOC")
sGroups <- list("FT", c("prPT", "rfPT", "arPT", "mrPT", "erPT"))
names(sGroups) <- sTypes
descrip <- "htseq_EdgeR_primary_HGSOC_vs_FT"
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC_CCNEamp_vs_HRD ###

#sTypes <- c("bothDrivers", "FT", "HRD", "CCNEamp", "unknown_driver")
#descrip <- "htseq_EdgeR_primary_HGSOC_CCNEamp_vs_HRD"
################################################################################

################################################################################
### htseq_EdgeR_HGSOC_drug_cats_vs_FT ###

#sTypes <- c("FT", "primary_resistant", "acquired_resistant", "drug_responders", 
#  "recurrent_ascites", "metastatic")
#sGroups <- list("FT", "prPT", "rfPT", "arPT", "mrPT", "erPT", "rcAF", "msST")
#names(sGroups) <- sTypes
#descrip <- "htseq_EdgeR_HGSOC_drug_cats_vs_FT"
################################################################################

################################################################################
### SalmonTE_primary_HGSOC_CCNEamp_vs_HRD ###

#sTypes <- c("bothDrivers", "FT", "HRD", "CCNEamp", "unknown_driver")
#descrip <- "SalmonTE_primary_HGSOC_CCNEamp_vs_HRD"
################################################################################

################################################################################
### htseq_EdgeR_primary_HGSOC_vs_FT_with_LINE1_silencers ###
#sTypes <- c("FT", "HGSOC")
#sGroups <- list("FT", c("prPT", "rfPT", "arPT", "mrPT", "erPT"))
#names(sGroups) <- sTypes
#descrip <- "htseq_EdgeR_primary_HGSOC_vs_FT_with_LINE1_silencers"
################################################################################

# define comparison parameters:
primaryOnly <- TRUE
cat_by_driver <- FALSE
EDAnormalise <- FALSE
count_tool <- "EdgeR"
customSamples <- FALSE

# define custom samples if needed:
#cus <- c("")

# specify what combination of repeat genes (repeats) and other genes,
# (all, both, other) should contribute to the results:
resultTypes <- c("repeats", "all", "both", "other")

# define sample group to use as control:
ctl <- "FT"

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.05
FCthresh <- 1

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
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                        expName, "/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                  "/plots/DEplots/", descrip, "/")

system(paste0("mkdir -p ", plotDir))
system(paste0("mkdir -p ", newRobjectDir))


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

# checkpoint to ensure Counts was loaded effectively:
if ( exists("Counts") ) {
    
  # if necessary, select primary samples only:
  if ( primaryOnly == TRUE) {
      Counts <- Counts[,grep("PT|FT", colnames(Counts))]
    }

  # select custom samples:
  if (customSamples) {
    Counts <- Counts[,colnames(Counts) %in% cus]
  }

  # re-categorize samples as HRD, CCNE_amp, both_drivers or unknown_drivers:
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
          colnames(Counts)[i] <- paste0(colnames(Counts)[i], "_unknown_drivers")
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
  
  if ( EDAnormalise == TRUE ) {
    # perform between lane full normalisation:
    nSet <- betweenLaneNormalization(set, which="full",offset=TRUE)
    
    # design matrix labelling all sample types:
    design <- model.matrix(~0+typeF, data=pData(nSet))
    
    # calculate the deviance residuals from a first-pass GLM regression of the counts he co-variates of interest (p8 RUVseq manual):
    # estimate dispersion:
    disp <- estimateGLMCommonDisp(counts(nSet),
                                  design, offset=-offst(nSet))
    
    # adjust values using dispersion:
    fit <- glmFit(counts(nSet), design, disp, offset=-offst(nSet))
    
    saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "DEfit.rds"))
    save.image(paste0(newRobjectDir, "/", Type, "DEdone.rds"))
    
  } else {
    
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
    
    # calculate normalisation factors and create post-norm RLE and PCA plots:
    if ( file.exists(paste0(RobjectDir, "/edgeRnorms.rds")) ) {
      norms <- readRDS(paste0(RobjectDir, "/edgeRnorms.rds"))
    } else {
      for (n in 1:nrow(Counts)) {
        print(n)
        if (n==1) {
          norms <- t(as.matrix(y$samples$norm.factors))
        } else {
          norms <- rbind(norms, norms[1,])
        }
      }
      
      saveRDS(norms, paste0(RobjectDir, "/edgeRnorms.rds"))
    }
    
    colnames(norms) <- colnames(Counts) 
    nSet <- newSeqExpressionSet(Counts, offset = norms, phenoData = data.frame(typeF, row.names=colnames(Counts)))
    
    # design matrix labelling all sample types:
    design <- model.matrix(~0+typeF)
    
    # estimate dispersion:
    disp <- estimateDisp(y, design=design)
    
    # adjust values using dispersion:
    fit <- glmFit(disp, design=design, robust=TRUE)
    
    saveRDS(fit, file=paste0(newRobjectDir, "/", Type, "DEfit.rds"))
    save.image(paste0(newRobjectDir, "/", Type, "DEdone.rds"))
    
  }
  
  
  # create post-norm RLE plot:
  if (file.exists(paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_RLElaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_RLElaneNormGC.pdf"))
    plotRLE(nSet, outline=FALSE, ylim=c(-4, 4))
    dev.off()
  }
  
  # create post-norm PCA:
  if (file.exists(paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf"))) {
    print(paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf already exists, no need to create"))
  } else {
    print(paste0("Creating ", plotDir, "/", Type, "_pcalaneNormGC.pdf"))
    pdf(file = paste0(plotDir, "/", Type, "_pcalaneNormGC.pdf"), height = 15, width = 15)
    plotPCA(nSet, cex=0.7)
    dev.off()
  }
  
  # determine which column has FT control:
  ctlInd <- grep(ctl, colnames(design))
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
  
  writeLines("\n")
  cont <- readline("Check above values - ok to continue? (y/n)")
  
  if (cont == "y") {
    
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
        
        
        ##############################################################################
        ### 4. Create DE data frames for repeats:
        ##############################################################################
        
        # define repeat and sig DE repeat dfs:
        repGenes <- allGenes[grep("ENS",  rownames(allGenes), invert = T),]
        print(repGenes)
        
        if ( is.na(FCthresh) ) {
          sigGenes <- filter(repGenes, FDR < FDRthresh)
          repGenes$threshold <- as.factor(repGenes$FDR < FDRthresh)
        } else if ( is.na(FDRthresh) ) {
          sigGenes <- repGenes[
            (repGenes$logFC > FCthresh)|(repGenes$logFC < -(FCthresh)), ]
          repGenes$threshold <- as.factor( 
            (repGenes$logFC > FCthresh)|(repGenes$logFC < -(FCthresh))
          )
        } else {
          sigGenes <- filter(
            repGenes, 
            (FDR < FDRthresh & logFC < -(FCthresh))|(
              FDR < FDRthresh & logFC > FCthresh
            )
          )
          repGenes$threshold <- as.factor(
            (repGenes$FDR < FDRthresh & repGenes$logFC < -(FCthresh))|(
              repGenes$FDR <  FDRthresh & repGenes$logFC > FCthresh
            )
            )
        }
        
        sig_rep <- subset(repGenes, threshold == T)
        
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
        
        if (length(FCthresh) == 0) {
          sigGenes <- filter(allGenes, FDR < FDRthresh)
          allGenes$threshold <- as.factor(allGenes$FDR < FDRthresh)
        } else {
          sigGenes <- filter(allGenes, (FDR < FDRthresh & logFC < -(FCthresh))|(FDR < FDRthresh & logFC > FCthresh))
          allGenes$threshold <- as.factor((allGenes$FDR < FDRthresh & allGenes$logFC < -(FCthresh))|(allGenes$FDR <  FDRthresh & allGenes$logFC > FCthresh))
        }
        
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
        ### 6. Create volcano plots:
        ##############################################################################
        
        if ("repeats|both" %in% resultTypes) {
        
        lab <- rbind(rbind(sig, posGenes), negGenes)
        repGenes <- rbind(rbind(repGenes, posGenes), negGenes)
        lab$genes <- rownames(lab)
        
        if ("both" %in% resultTypes) {
        }
        
        }
        
        if ("all" %in% resultTypes) {
        }

        

        
        
        
        }
    } 
        
        
        
    
}

