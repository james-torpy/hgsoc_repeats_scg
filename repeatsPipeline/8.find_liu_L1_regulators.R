### 19.find_L1_regulators.R ###

# Finds suppressors from Liu et. al. 2018 that are significantly 
# upregulated when L1s are downregulated (FDR < 0.1) or 
# significantly downregulated when L1s are upregulated 
# (FDR < 0.1) in at least one DE comparison, and follow this 
# trend for all other DE comparisons for FDR < 0.3. Also must 
# follow what they are assigned as in Liu et. al. 2018. Also finds
# activators

library(tibble)
library(dplyr)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(edgeR)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
Type <- "custom3"
descrip <- "find_liu_L1_regulators"

################################################################################
### Options ###
################################################################################

################################################################################
### primary_HGSOC_vs_FT ###
sTypes <- list(c("FT", "HGSOC"))
sGroups <- list(list("FT", c("prPT", "rfPT", "arPT", "mrPT", "erPT")))
names(sGroups[[1]]) <- sTypes[[1]]
ctl <- list("FT")

### drug_cats_vs_primary_resistance ###
sTypes[[2]] <- c("primary_resistant", "acquired_resistant", "multiple_responders",
            "extreme_responders")
sGroups[[2]] <- list(c("prPT", "rfPT"), "arPT", "mrPT", "erPT")
names(sGroups[[2]]) <- sTypes[[2]]
ctl[[2]] <- list("primary_resistant")

### drug_cats_vs_acquired_resistance ###
sTypes[[3]] <- c("primary_resistant", "acquired_resistant", "multiple_responders",
             "extreme_responders")
sGroups[[3]] <- list(c("prPT", "rfPT"), "arPT", "mrPT", "erPT")
names(sGroups[[3]]) <- sTypes[[3]]
ctl[[3]] <- list("acquired_resistant")

names(sTypes) <- c("primary_HGSOC_vs_FT", "drug_cats_vs_primary_resistance",
                   "drug_cats_vs_acquired_resistance")
names(sGroups) <- c("primary_HGSOC_vs_FT", "drug_cats_vs_primary_resistance",
                   "drug_cats_vs_acquired_resistance")
names(ctl) <- c("primary_HGSOC_vs_FT", "drug_cats_vs_primary_resistance",
                   "drug_cats_vs_acquired_resistance")

################################################################################

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
refDir <- paste0(projectDir, "/RNA-seq/refs/")
rawDir <- paste0(projectDir, 
                 "/RNA-seq/raw_files/fullsamples/bowtell_primary/")
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/")
newRobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", descrip, "/")
tableDir <- paste0(resultsDir, "/R/tables/DEplots/", descrip, "/")

system(paste0("mkdir -p ", tableDir))
system(paste0("mkdir -p ", newRobjectDir))

# specify other genes to include if necessary:
other_df <- read.csv(paste0(refDir, "/liu2018_L1_regulators.csv"), header=T)[,1:3]

# reduce list of other genes:
# red <- c("RNASEH2A", "RAD54L", "RAD51", "UBE2T", "MCM2", "MCM4", "MCM3", "MCM6", 
#          "MCM8", "CWC15", "FANCE", "XRCC2", "FANCB", "UQCRH", "ATP5L", "BRCA2", 
#          "FANCI", "RAB36", "HECTD1", "LIG4", "DLG4", "EHMT1", "CHST14", "USP15", 
#          "TSPYL1", "ALKBH1", "XPNPEP1")
# other_df <- other_df[other_df$symbol %in% red,]

colnames(other_df) <- c("ensembl_id", "symbol", "type")
otherIDs <- other_df$ensembl_id
otherSym <- other_df$symbol
suppressor_df <- other_df[other_df$type == "suppressor",]
activator_df <- other_df[other_df$type == "activator",]


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

  
if ( !file.exists(paste0(RobjectDir, "/", Type, "_EdgeR",
                         "_counts.rds")) ) {
  
  writeLines("\n")
  print("EdgeR counts data frame does not exist, creating now...")
  
  custom3Counts <- readRDS(paste0(RobjectDir, "/", Type, 
                                  "_allcounts.htseq.rds"))
  gcCounts <- readRDS(paste0(RobjectDir, "/gc_allcounts.htseq.rds"))
  
  # append gcCounts to custom3Counts:
  Counts <- counts_bind(custom3Counts, gcCounts)
  
  saveRDS(Counts, paste0(RobjectDir, "/", Type, "_EdgeR", 
                         "_counts.rds"))
} else {
  
  print("Loading EdgeR counts data frame...")
  Counts <- readRDS(paste0(RobjectDir, "/", Type, "_EdgeR", 
                           "_counts.rds"))
}

# select primary samples only:
Counts <- Counts[,grep("PT|FT", colnames(Counts))]

# eliminate lowly expressed genes (rows where there are less than 3 counts 
# where df > 4):
print(paste0("No. rows before filtering is: ", nrow(Counts)))
Counts <- Counts %>%
  rownames_to_column('gene_id') %>%
  dplyr::filter(rowSums(Counts > 5) >= (ncol(Counts)/3)) %>%
  column_to_rownames('gene_id')
print(paste0("No. rows after  filtering: ", nrow(Counts)))


################################################################################
### 2. Group counts ###
################################################################################
counts_temp <- Counts

j=1
DE_regulator_list <- lapply(sGroups, function(x) {
  
  # change sample names according to grouping:
  for (i in 1:length(x)) {
    for (n in x[[i]]) {
      colnames(counts_temp) <- gsub(n, names(x)[i], colnames(counts_temp))
    }
  }
  
  # remove samples not belonging to any group:
  counts_temp <- counts_temp[,gsub("AOCS_[0-9][0-9][0-9]_", "", colnames(counts_temp)) %in% names(x)]
  
  # change the order of columns of counts_temp to alphabetical order:
  counts_temp <- counts_temp[,order(
    gsub(
      "AOCS.*_[0-9][0-9][0-9]_", "", colnames(counts_temp)
    )
  )]
  
  # define sample groups:
  splt <- unlist(
    lapply(
      split(
        colnames(counts_temp), gsub(
          "AOCS.*_[0-9][0-9][0-9]_", "", colnames(counts_temp)
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
  levels(typeF) <- sTypes[[j]]
  
  sampleNos <- unlist(
    lapply(
      split(
        colnames(counts_temp), gsub(
          "\\.1", "", gsub(
            "AOCS.*_[0-9][0-9][0-9]_", "", colnames(counts_temp)
          )
        )
      ), length
    )
  )
  
  # convert Counts1 into SeqExpressionSet - elements need to be delisted and 
  # changed to integers first:
  counts_temp <- apply(counts_temp, 2, unlist)
  storage.mode(counts_temp) <- "integer"
  
  
  ##############################################################################
  ### 3. Perform normalisation and DE:
  ##############################################################################
  
  # perform between lane full normalisation:
  y <- DGEList(counts = counts_temp, group = typeF)
  
  # normalise for library size:
  y <- calcNormFactors(y)
  
  # design matrix labelling all sample types:
  design <- model.matrix(~0+typeF)
  colnames(design) <- gsub("typeF", "", colnames(design))
  
  # estimate dispersion:
  disp <- estimateDisp(y, design=design)
  
  # adjust values using dispersion:
  fit <- glmFit(disp, design=design, robust=TRUE)
  
  # determine which column has control:
  ctlInd <- which(colnames(design)==ctl[[j]])
  
  con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  
  # put sTypes in alphabetical order:
  sTypes[[j]] <- sTypes[[j]][order(sTypes[[j]])]
  
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
      
      comp <- paste0(sTypes[[j]][i], "_vs_", ctl[[j]])
      print(paste0("Finding regulators from ", comp))
      
      # perform likelihood ratio test:
      con[i] <- 1
      lrt <- glmLRT(fit, contrast = con)
      
      # determine the top DE genes:
      topTags(lrt)
      
      # fetch summary of differentially expressed genes (those  with FDR =< 0.05:
      DEs <- summary(result <- decideTestsDGE((lrt)))
      
      # fetch all gene DE info, 
      allGenes <- as.data.frame(topTags(lrt, n=Inf))
      
      
      ##############################################################################
      ### 4. Determine whether L1s are down, up or not differentially regulated:
      ##############################################################################
      
      # isolate L1s with FDR < 0.1:
      L1s <- allGenes[grep("^L1", rownames(allGenes)),]
      L1s <- L1s[L1s$FDR < 0.1,]
      
      # determine direction of L1 differential expression if there is DE:
      if ( nrow(L1s) > 0 ) {
        if ( mean(L1s$logFC) > 0 ) {
          L1_dir <- "up"
        } else if ( mean(L1s$logFC) < 0 ) {
          L1_dir <- "down"
        } else {
          L1_dir <- "non_DE"
        }
      } else {
        L1_dir <- "non_DE"
      }
      
      ##############################################################################
      ### 5. Create DE data frames for gencode genes:
      ##############################################################################
      
      if ( L1_dir != "non_DE" ) {
        
        # annotate allGenes with entrez ids and symbols in separate columns:
        egENSEMBL <- toTable(org.Hs.egENSEMBL)
        egSYMBOL <- toTable(org.Hs.egSYMBOL)
        
        allGenes$gene_id <- egENSEMBL$gene_id[match(rownames(allGenes), egENSEMBL$ensembl_id)]
        allGenes$symbol <- egSYMBOL$symbol[match(allGenes$gene_id, egSYMBOL$gene_id)]
        
        # isolate Liu 2018 L1 regulators:
        liu <- allGenes[rownames(allGenes) %in% otherIDs,] 
        
        # create suppressor and activator lists, fdr < 0.3:
        trend_liu <- liu[liu$FDR < 0.3,]
        
        if ( L1_dir == "down" ) {
          
          if ( !exists("suppressors_fdr0.3") ) {
            
            suppressors_fdr0.3 <- trend_liu[trend_liu$logFC > 0,]
            suppressors_fdr0.3 <- suppressors_fdr0.3[rownames(suppressors_fdr0.3) %in% suppressor_df$ensembl_id,]
            activators_fdr0.3 <- trend_liu[trend_liu$logFC < 0,]
            activators_fdr0.3 <- activators_fdr0.3[rownames(activators_fdr0.3) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(suppressors_fdr0.3) > 0 ) {
              suppressors_fdr0.3$comp <- comp
              suppressors_fdr0.3$L1_dir <- L1_dir
            }
            if ( nrow(activators_fdr0.3) > 0 ) {
              activators_fdr0.3$comp <- comp
              activators_fdr0.3$L1_dir <- L1_dir
            }
            
          } else {
            
            sup0.3 <- trend_liu[trend_liu$logFC > 0,]
            sup0.3[rownames(sup0.3) %in% suppressor_df$ensembl_id,]
            act0.3 <- trend_liu[trend_liu$logFC < 0,]
            act0.3 <- act0.3[rownames(act0.3) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(sup0.3) > 0 ) {
              sup0.3$comp <- comp
              sup0.3$L1_dir <- L1_dir
            }
            if ( nrow(act0.3) > 0 ) {
              act0.3$comp <- comp
              act0.3$L1_dir <- L1_dir
            }
            
            # add new regulators to existing dfs:
            suppressors_fdr0.3 <- rbind(suppressors_fdr0.3, sup0.3)
            activators_fdr0.3 <- rbind(activators_fdr0.3, act0.3)
            
          }
          
          # create suppressor and activator lists, fdr < 0.1:
          trend_liu <- liu[liu$FDR < 0.1,]
          if ( !exists("suppressors_fdr0.1") ) {
            
            suppressors_fdr0.1 <- trend_liu[trend_liu$logFC > 0,]
            suppressors_fdr0.1 <- suppressors_fdr0.1[rownames(suppressors_fdr0.1) %in% suppressor_df$ensembl_id,]
            activators_fdr0.1 <- trend_liu[trend_liu$logFC < 0,]
            activators_fdr0.1 <- activators_fdr0.1[rownames(activators_fdr0.1) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(suppressors_fdr0.1) > 0 ) {
              suppressors_fdr0.1$comp <- comp
              suppressors_fdr0.1$L1_dir <- L1_dir
            }
            if ( nrow(activators_fdr0.1) > 0 ) {
              activators_fdr0.1$comp <- comp
              activators_fdr0.1$L1_dir <- L1_dir
            }
            
          } else {
            
            sup0.1 <- trend_liu[trend_liu$logFC > 0,]
            sup0.1[rownames(sup0.1) %in% suppressor_df$ensembl_id,]
            act0.1 <- trend_liu[trend_liu$logFC < 0,]
            act0.1 <- act0.1[rownames(act0.1) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(sup0.1) > 0 ) {
              sup0.1$comp <- comp
              sup0.1$L1_dir <- L1_dir
            }
            if ( nrow(act0.1) > 0 ) {
              act0.1$comp <- comp
              act0.1$L1_dir <- L1_dir
            }
            
            # add new regulators to existing dfs:
            suppressors_fdr0.1 <- rbind(suppressors_fdr0.1, sup0.1)
            activators_fdr0.1 <- rbind(activators_fdr0.1, act0.1)
            
          }
          
        } else if ( L1_dir == "up" ) {
          
          if ( !exists("suppressors_fdr0.3") ) {
            
            suppressors_fdr0.3 <- trend_liu[trend_liu$logFC < 0,]
            suppressors_fdr0.3 <- suppressors_fdr0.3[rownames(suppressors_fdr0.3) %in% suppressor_df$ensembl_id,]
            activators_fdr0.3 <- trend_liu[trend_liu$logFC > 0,]
            activators_fdr0.3 <- activators_fdr0.3[rownames(activators_fdr0.3) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(suppressors_fdr0.3) > 0 ) {
              suppressors_fdr0.3$comp <- comp
              suppressors_fdr0.3$L1_dir <- L1_dir
            }
            if ( nrow(activators_fdr0.3) > 0 ) {
              activators_fdr0.3$comp <- comp
              activators_fdr0.3$L1_dir <- L1_dir
            }
            
          } else {
            
            sup0.3 <- trend_liu[trend_liu$logFC < 0,]
            sup0.3 <- sup0.3[rownames(sup0.3) %in% suppressor_df$ensembl_id,]
            act0.3 <- trend_liu[trend_liu$logFC > 0,]
            act0.3 <- act0.3[rownames(act0.3) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(sup0.3) > 0 ) {
              sup0.3$comp <- comp
              sup0.3$L1_dir <- L1_dir
            }
            if ( nrow(act0.3) > 0 ) {
              act0.3$comp <- comp
              act0.3$L1_dir <- L1_dir
            }
            
            # add new regulators to existing dfs:
            suppressors_fdr0.3 <- rbind(suppressors_fdr0.3, sup0.3)
            activators_fdr0.3 <- rbind(activators_fdr0.3, act0.3)
            
          }
          
          # create suppressor and activator lists, fdr < 0.1:
          trend_liu <- liu[liu$FDR < 0.1,]
          if ( !exists("suppressors_fdr0.1") ) {
            
            suppressors_fdr0.1 <- trend_liu[trend_liu$logFC < 0,]
            suppressors_fdr0.1 <- suppressors_fdr0.1[rownames(suppressors_fdr0.1) %in% suppressor_df$ensembl_id,]
            activators_fdr0.1 <- trend_liu[trend_liu$logFC > 0,]
            activators_fdr0.1 <- activators_fdr0.1[rownames(activators_fdr0.1) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(suppressors_fdr0.1) > 0 ) {
              suppressors_fdr0.1$comp <- comp
              suppressors_fdr0.1$L1_dir <- L1_dir
            }
            if ( nrow(activators_fdr0.1) > 0 ) {
              activators_fdr0.1$comp <- comp
              activators_fdr0.1$L1_dir <- L1_dir
            }
            
          } else {
            
            sup0.1 <- trend_liu[trend_liu$logFC < 0,]
            sup0.1[rownames(sup0.1) %in% suppressor_df$ensembl_id,]
            act0.1 <- trend_liu[trend_liu$logFC > 0,]
            act0.1 <- act0.1[rownames(act0.1) %in% activator_df$ensembl_id,]
            
            # add comparison info and L1 direction columns to dfs:
            if ( nrow(sup0.1) > 0 ) {
              sup0.1$comp <- comp
              sup0.1$L1_dir <- L1_dir
            }
            if ( nrow(act0.1) > 0 ) {
              act0.1$comp <- comp
              act0.1$L1_dir <- L1_dir
            }
            
            # add new regulators to existing dfs:
            suppressors_fdr0.1 <- rbind(suppressors_fdr0.1, sup0.1)
            activators_fdr0.1 <- rbind(activators_fdr0.1, act0.1)
            
          }
        }
      }
    }
    con <- c(rep(0, (ctlInd - 1) ), -1, rep(0, (ncol(design) - ctlInd)))
  }
  
  res <- list(suppressors_fdr0.3, activators_fdr0.3, suppressors_fdr0.1, activators_fdr0.1)
  names(res) <- c("suppressors_fdr0.3", "activators_fdr0.3", "suppressors_fdr0.1", "activators_fdr0.1")
  
  j <<- j+1
  return(res)
  
})

for ( m in 1:length(DE_regulator_list) ) {
  if (m==1) {
    sup0.3_df <- DE_regulator_list[[m]][[1]]
    act0.3_df <- DE_regulator_list[[m]][[2]]
    sup0.1_df <- DE_regulator_list[[m]][[3]]
    act0.1_df <- DE_regulator_list[[m]][[4]]
  } else {
    sup0.3_df <- rbind(sup0.3_df, DE_regulator_list[[m]][[1]])
    act0.3_df <- rbind(act0.3_df, DE_regulator_list[[m]][[2]])
    sup0.1_df <- rbind(sup0.1_df, DE_regulator_list[[m]][[3]])
    act0.1_df <- rbind(act0.1_df, DE_regulator_list[[m]][[4]])
  }
}
final_sup <- sup0.1_df[rownames(sup0.1_df) %in% suppressor_df$ensembl_id,]
final_sup <- as.data.frame(unique(final_sup$symbol[!(final_sup$symbol %in% act0.3_df$symbol)]), 
                           row.names = row.names(final_sup))
colnames(final_sup) <- "symbol"

final_act <- act0.1_df[rownames(act0.1_df) %in% activator_df$ensembl_id,]
final_act <- as.data.frame(unique(final_act$symbol[!(final_act$symbol %in% sup0.3_df$symbol)]), 
                           row.names = row.names(final_act))
colnames(final_act) <- "symbol"

# sanity checks:
# find if there are genes in final_sup not in suppressor_df:
discrep_sup <- final_sup$symbol[!(final_sup$symbol %in% suppressor_df$symbol)]

# for discrepencies check if the ensembl IDs match:
d_ids <- rownames(sup0.1_df)[sup0.1_df$symbol %in% discrep_sup]
d_ids %in% suppressor_df$ensembl_id

# do above for activators:
discrep_act <- final_act$symbol[!(final_act$symbol %in% activator_df$symbol)]
d_ids <- rownames(act0.1_df)[act0.1_df$symbol %in% discrep_act]
d_ids %in% activator_df$ensembl_id

# save as tables:
final_sup$ensembl_id <- rownames(final_sup)
write.table(final_sup, paste0(tableDir, "/verified_liu_L1_suppressors.txt"), 
                sep = "\t", row.names = F, col.names = T, quote = F)
final_act$ensembl_id <- rownames(final_act)
write.table(final_act, paste0(tableDir, "/verified_liu_L1_activators.txt"), 
            sep = "\t", row.names = F, col.names = T, quote = F)

