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
 descrip <- "htseq_EdgeR_primary_HGSOC_vs_FT_with_LINE1_repressors"
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
resultTypes <- c("repeats", "gc", "both", "other")

# define sample group to use as control:
ctl <- "FT"

# specify what FDR and log2 fold change thresholds to use:
FDRthresh <- 0.1
FCthresh <- 1

# specify control genes to include:
posGeneIDs <- c("ENSG00000111640", "ENSG00000196776")
posGeneNames <- c("GAPDH", "CD47")
negGeneIDs <- c("ENSG00000075624", "ENSG00000169919")
negGeneNames <- c("beta-actin", "GUSB")

# specify other genes to include if necessary:
otherIDs <- c(# chromatin remodellers:
            "ENSG00000128908", "ENSG00000183495", "ENSG00000080603", "ENSG00000153922", 
            "ENSG00000173575", "ENSG00000170004", "ENSG00000111642", "ENSG00000124177", 
            "ENSG00000171316", "ENSG00000100888", "ENSG00000177200", "ENSG00000116254", 
            "ENSG00000080503", "ENSG00000127616", "ENSG00000153147", "ENSG00000102038",
            "ENSG00000076108")

otherSym <- c (# chromatin remodellers:
              "INO80", "EP400", "SRCAP", "CHD1", 
              "CHD2", "CHD3", "CHD4", "CHD6", 
              "CHD7", "CHD8", "CHD9", "CHD5", 
              "SMARCA2", "SMARCA4", "SMARCA5", "SMARCA1",
              "BAZ2A")

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
    save.image(paste0(newRobjectDir, "/", Type, "DEfit_done.rds"))
    
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
        
        ######
        colnames(gcCounts)[1] <- "ensembl_id"
        gcCounts$gene_id <- egENSEMBL$gene_id[match(gcCounts$ensembl_id, egENSEMBL$ensembl_id)]
        gcCounts$symbol <- egSYMBOL$symbol[match(gcCounts$gene_id, egSYMBOL$gene_id)]
        outs <- gcCounts[gcCounts$symbol %in% rownames(outliers),]
        apply(subset(outs, select=-c(ensembl_id, gene_id, symbol)), 1, hist)
        low_counts <- outliers[c("ANKRD19P", "BEST1", "CXCR5", "DND1", "FOXH1", "TUBB2B"),]
        lows <- gcCounts[gcCounts$symbol %in% rownames(low_counts),]
        ######

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
        ctlGenes$genes<- rownames(ctlGenes)


        ##############################################################################
        ### 7. Create volcano plots:
        ##############################################################################
        
        if ("both" %in% resultTypes) {
        
          lab <- rbind(sig_rep, ctlGenes)
          lab$genes <- rownames(lab)
          bothGenes <- rbind(rbind(repGenes, gcGenes), ctlGenes)

          # combine 'threshold' and 'type' columns:
          bothGenes$type_thresh <- paste0(bothGenes$type, "_", bothGenes$threshold)
          lab$type_thresh <- paste0(lab$type, "_", lab$threshold)

          # plot on volcano plot:
          p <- ggplot(data=bothGenes, aes( x=logFC, y=-log10(FDR), color=type_thresh))
          p <- p + geom_point(data=bothGenes)
          p <- p + geom_text_repel(data=lab, aes(label=genes))
          p <- p + theme(legend.position =  "none")
          p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
          # key for colours = c("neg_ctls", "pos_ctls", "neg_reps", "gc", "gc", "pos_reps")
          #p <- p + scale_colour_manual(values = c("#114477", "firebrick4", 
#            "dodgerblue1", "#C6C6C6", "#C6C6C6", "firebrick1"))
          p <- p + scale_colour_manual(values = c("darkorchid1", "darkorchid4", 
            "darkolivegreen3", "#C6C6C6", "#C6C6C6", "forestgreen"))
          #p <- p +  xlim(c(-2, 2))
          if (length(FCthresh) == 0) {
            if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_all_genes.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_all_genes.pdf"))
              p
            } else {
              print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_all_genes.pdf"))
              pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, "_all_genes.pdf"))
              print(p)
              dev.off()
            }
          } else {
            if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_all_genes.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_all_genes.pdf already exists"))
              p
            } else {
              print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, "_all_genes.pdf"))
              pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_all_genes.pdf"))
              print(p)
              dev.off()
            }
          
        }
        
        if ("repeats" %in% resultTypes) {
          lab <- rbind(sig_rep, ctlGenes)
          lab$genes <- rownames(lab)
          repGenes <- rbind(repGenes, ctlGenes)

          # combine 'threshold' and 'type' columns:
          repGenes$type_thresh <- paste0(repGenes$type, "_", repGenes$threshold)
          lab$type_thresh <- paste0(lab$type, "_", lab$threshold)

          # plot on volcano plot:
          p <- ggplot(data=repGenes, aes( x=logFC, y=-log10(FDR), color=type_thresh))
          p <- p + geom_point(data=repGenes)
          p <- p + geom_text_repel(data=lab, aes(label=genes))
          p <- p + theme(legend.position =  "none")
          p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
          # key for colours = c("neg_ctls", "pos_ctls", "neg_gc", "pos_gc")
#          p <- p + scale_colour_manual(values = c("#114477", "firebrick4", 
#            "dodgerblue1", "firebrick1"))
          p <- p + scale_colour_manual(values = c("darkorchid1", "darkorchid4", 
            "darkolivegreen3", "#C6C6C6", "#C6C6C6", "forestgreen"))
          #p <- p +  xlim(c(-2, 2))
          if (length(FCthresh) == 0) {
            if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_reps.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_reps.pdf"))
              p
            } else {
              print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_reps.pdf"))
              pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, "_reps.pdf"))
              print(p)
              dev.off()
            }
          } else {
            if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf already exists"))
              p
            } else {
              print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf"))
              pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_reps.pdf"))
              print(p)
              dev.off()
            }
          } 

        } else if ("gc" %in% resultTypes) {

          lab <- rbind(sig_gc, ctlGenes)
          lab$genes <- rownames(lab)
          bothGenes <- rbind(gcGenes, ctlGenes)

          # combine 'threshold' and 'type' columns:
          gcGenes$type_thresh <- paste0(gcGenes$type, "_", gcGenes$threshold)
          lab$type_thresh <- paste0(lab$type, "_", lab$threshold)

          # plot on volcano plot:
          p <- ggplot(data=gcGenes, aes( x=logFC, y=-log10(FDR), color=type_thresh))
          p <- p + geom_point(data=gcGenes)
          p <- p + geom_text_repel(data=lab, aes(label=genes))
          p <- p + theme(legend.position =  "none")
          # key for colours = c("neg_ctls", "pos_ctls", "neg_gc", "pos_gc")
#          p <- p + scale_colour_manual(values = c("#114477", "firebrick4", 
#            "dodgerblue1", "firebrick1"))
          p <- p + scale_colour_manual(values = c("darkorchid1", "darkorchid4", 
            "darkolivegreen3", "#C6C6C6", "#C6C6C6", "forestgreen"))
          #p <- p +  xlim(c(-2, 2))
          if (length(FCthresh) == 0) {
            if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_gc_genes.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_gc_genes.pdf"))
              p
            } else {
              print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_gc_genes.pdf"))
              pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, "_gc_genes.pdf"))
              print(p)
              dev.off()
            }
          } else {
            if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_gc_genes.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_gc_genes.pdf already exists"))
              p
            } else {
              print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, "_gc_genes.pdf"))
              pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_gc_genes.pdf"))
              print(p)
              dev.off()
            }
          }
        } else if ("other" %in% resultTypes) {

          otherGenes <- allGenes[otherIDs,]
          otherGenes$type <- "other"
          otherGenes$genes <- rownames(otherGenes)
          
          sig_other <- sig_gc[sig_gc$symbol %in% otherSym,]
          sig_other$type <- "other"
          sig_other$genes <- rownames(sig_other)
          lab <- rbind(sig_other, ctlGenes)
          lab$genes <- rownames(lab)
          
          otherGenes <- rbind(otherGenes, ctlGenes)

          # combine 'threshold' and 'type' columns:
          otherGenes$type_thresh <- paste0(otherGenes$type, "_", otherGenes$threshold)
          lab$type_thresh <- paste0(lab$type, "_", lab$threshold)

          # plot on volcano plot:
          p <- ggplot(data=otherGenes, aes( x=logFC, y=-log10(FDR), color=type_thresh))
          p <- p + geom_point(data=otherGenes)
          p <- p + geom_text_repel(data=lab, aes(label=symbol))
          p <- p + theme(legend.position =  "none")
          p <- p + labs(x="log2 fold change vs FT control", y="-log10 FDR")
          # key for colours = c("neg_ctls", "pos_ctls", "neg_gc", "pos_gc")
          #p <- p + scale_colour_manual(values = c("#114477", "firebrick4", 
           # "dodgerblue1", "firebrick1", "forestgreen"))
          #p <- p + scale_colour_manual(values = c("darkorchid1", "darkorchid4", 
          #  "darkolivegreen3", "#C6C6C6", "#C6C6C6", "forestgreen"))
          #p <- p +  xlim(c(-2, 2))
          if (length(FCthresh) == 0) {
            if (file.exists(paste0(plotDir,   "/", Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_other_genes.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR_",   FDRthresh, "_", comp, "_other_genes.pdf"))
              p
            } else {
              print(paste0("Creating  ",plotDir, "/", Type,    "_volcano_FDR_", FDRthresh, "_", comp, "_other_genes.pdf"))
              pdf(file = paste0(plotDir, "/",   Type,  "_volcano_FDR_",  FDRthresh, "_", comp, "_other_genes.pdf"))
              print(p)
              dev.off()
            }
          } else {
            if (file.exists(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_other_genes.pdf"))) {
              print(paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_other_genes.pdf already exists"))
              p
            } else {
              print(paste0("Creating  ", plotDir, "/",  Type,  "_volcano_FDR", FDRthresh, "_FC", FCthresh, "_", comp, "_other_genes.pdf"))
              pdf(file = paste0(plotDir, "/",  Type,  "_volcano_FDR",   FDRthresh, "_FC", FCthresh, "_", comp, "_other_genes.pdf"))
              print(p)
              dev.off()
            }
          }
        }
        con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
      }
    } 

    ##############################################################################
    ### 8. Save all results dfs ###
    ##############################################################################
    
    # remove the NULL list dfs created when avoiding clt vs ctl:
    if (length(sTypes)>2) {
      allReps <- allReps[-ctlInd]
      sigReps <- sigReps[-ctlInd]
      # name the list elements:
      names(allReps) <- paste0(sTypes[-ctlInd], "_vs_", ctl)
      names(sigReps) <- paste0(sTypes[-ctlInd], "_vs_", ctl)
    }
    
    if (file.exists(paste0(newRobjectDir, "/", Type, "_DEallReps.rds"))) {
      print(paste0(newRobjectDir, "/", Type, "_DEallReps.rds already exists"))
    } else {
      saveRDS(allReps, file=paste0(newRobjectDir, "/", Type, "_DEallReps.rds"))
    }
    
    if (file.exists(paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))) {
      print(paste0(newRobjectDir, "/", Type, "_DEsigReps.rds already exists"))
    } else {
      saveRDS(sigReps, file=paste0(newRobjectDir, "/", Type, "_DEsigReps.rds"))
    }
    
    if (file.exists(paste0(newRobjectDir, "/", Type, "_DEallGenes.rds"))) {
      print(paste0(newRobjectDir, "/", Type, "_DEallGenes.rds already exists"))
    } else {
      saveRDS(allGenes, file=paste0(newRobjectDir, "/", Type, "_DEallGenes.rds"))
    }
    
    if (file.exists(paste0(newRobjectDir, "/", Type, "_DEsigGenes.rds"))) {
      print(paste0(newRobjectDir, "/", Type, "_DEsigGenes.rds already exists"))
    } else {
      saveRDS(sigGenes, file=paste0(newRobjectDir, "/", Type, "_DEsigGenes.rds"))
    }   

  # if values of con, design are not correct, print error message: 
  } else {
    print("Check values and run again")
  }

# if Counts was not loaded, print error message:
} else {
  print("Error - Counts did not load")     
}
