### 20.find_correlated_L1_regulators.R ###

# Takes verified regulators from Liu et. al. 2018 and tests
# for correlations between CPM of regulator and mean differentially
# expressed L1 CPMR:

library(tibble)
library(dplyr)
library(Rmisc)
library(org.Hs.eg.db)

# define starting variables:
project <- "hgsoc_repeats"
expName <- "exp9"
Type <- "custom3"
descrip <- "find_liu_L1_regulators"

# define directories:
homeDir <- "/Users/jamestorpy/clusterHome/"
#homeDir <- "/share/ScratchGeneral/jamtor/"
projectDir <- paste0(homeDir, "/projects/", project)
resultsDir <- paste0(projectDir, "/RNA-seq/results")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/",
                     expName, "/")
refDir <- paste0(projectDir, "/RNA-seq/refs/")
tableDir <- paste0(resultsDir, "/R/", expName,
                   "/tables/DEplots/", descrip, "/")
plotDir <- paste0(resultsDir, "/R/", expName,
                   "/plots/DEplots/", descrip, "/") 

system(paste0("mkdir -p ", tableDir))
system(paste0("mkdir -p ", plotDir))


################################################################################
### 1. Calculate CPMs and isolate verified L1 regulators and differentially 
# expressed L1s ###
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

# calculate total repeat count size:
cSizes <- apply(Counts, 2, sum)

# calculate CPMs
CPM <- as.data.frame(t(t(Counts)/cSizes)*1000000)
#CPM <- CPM[,grep("FT", colnames(CPM), invert = T)]

# load L1 regulators:
sups <- read.table(paste0(tableDir, "/verified_liu_L1_suppressors.txt"),
                   sep = "\t", header = T)
acts <- read.table(paste0(tableDir, "/verified_liu_L1_activators.txt"),
                                sep = "\t", header = T)

# reg_df <- read.csv(paste0(refDir, "/liu2018_L1_regulators.csv"), header=T)[,1:3]
# sups <- reg_df[reg_df$effect == "suppressor",]
# acts <- reg_df[reg_df$effect == "activator",]

# isolate suppressor and activator CPMs:
sup_CPM <- CPM[as.character(sups$ensembl_id),]
act_CPM <- CPM[as.character(acts$ensembl_id),]

# load differentially expressed L1s:
DE <- readRDS(file = paste0(RobjectDir, "/DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/custom3_DEsigReps.rds"))
repeatList <- rownames(DE[[1]])
DE_L1s <- grep("^L1", repeatList, value = T)

# isolate CPMs for differentially expressed L1s:
L1_CPM <- CPM[DE_L1s,]

# calculate mean L1 value for each sample:
L1_CPM_mean <- apply(L1_CPM, 2, mean)

# assess distribution of mean L1 CPM:
histogram(L1_CPM_mean)


################################################################################
### 2. For each suppressor and activator, calculate correlation with 
# mean CPM of DE L1s and save scatterplots ###
################################################################################

# build ID and symbol data frames:
egENSEMBL <- toTable(org.Hs.egENSEMBL)
egSYMBOL <- toTable(org.Hs.egSYMBOL)

for ( i in 1:nrow(sup_CPM) ) {
  print(paste0("Finding correlation for ", rownames(sup_CPM[i,]),
               " and mean_L1_CPM"))
  
  # create data frame with suppressor CPMs and mean L1 CPMs:
  df <- data.frame(as.numeric(sup_CPM[i,]), as.numeric(L1_CPM_mean))
  colnames(df) <- c("suppressor_CPM", "mean_L1_CPM")
  df <- df[apply(df, 1, function(row) all(row !=0 )),]
  
  # calculate Spearman correlation:
  sup_cor <- cor(df, method = "spearman")[2,1]
  
  if (i==1) {
    sup_cors <- list(list(df, sup_cor))
    names(sup_cors)[i] <- rownames(sup_CPM[i,])
  } else {
    sup_cors[[i]] <- list(df, sup_cor)
    names(sup_cors)[i] <- rownames(sup_CPM[i,])
  }
  
  # determine symbol corresponding to ID:
  m <- match(rownames(sup_CPM[i,]), egENSEMBL$ensembl_id)
  gene_id <- egENSEMBL$gene_id[m]
  m <- match(gene_id, egSYMBOL$gene_id)
  sym <- egSYMBOL$symbol[m]
  
  # if correlation < -0.5, save as scatterplot:
  if ( sup_cor < -0.3 ) {
    log_df <- log10(df)
    
    reg <- lm(mean_L1_CPM~suppressor_CPM, data = log_df)
    
    # reg <- lm(suppressor_CPM~mean_L1_CPM, data = log_df)
    xmax <- max(log_df$suppressor_CPM)
    ymax <- max(log_df$mean_L1_CPM)
    
    pdf(paste0(plotDir, "/", sym, 
               "_suppressor_correlation_with_mean_L1_CPM_exp.pdf"))
    par(mar = c(5,5,3,3))
    plot(log_df, ann = F)
    abline(reg, col="red")
    text(x=xmax-(0.1), y=ymax-(0.1), paste0("R=", round(sup_cor, 3)))
    title(main=paste0("mean L1 vs ", sym), xlab = paste0(sym, 
      " suppressor CPM"), ylab = paste0("mean L1 CPM"))
    dev.off()
    
    if ( !(exists("correlated_sups_0.5")) ) {
      correlated_sups_0.5 <- c(paste0("mean_L1_vs_", sym))
    } else {
      correlated_sups_0.5 <- c(correlated_sups_0.5, paste0("mean_L1_vs_", sym))
    }
  }
  
  # if correlation < -0.3, save as table:
  if ( sup_cor < -0.3 ) {
    if ( !(exists("correlated_sups_0.3")) ) {
      correlated_sups_0.3 <- c(paste0("mean_L1_vs_", sym))
    } else {
      correlated_sups_0.3 <- c(correlated_sups_0.3, paste0("mean_L1_vs_", sym))
    }
  }
}
write.table(correlated_sups_0.3, paste0(tableDir, 
                                        "/L1_suppressors_correlated_with_mean_L1_CPMs_0.3.txt"), quote = F, sep = "\n")
write.table(correlated_sups_0.5, paste0(tableDir, 
                                        "/L1_suppressors_correlated_with_mean_L1_CPMs_0.5.txt"), quote = F, sep = "\n")


for ( i in 1:nrow(act_CPM) ) {
  print(paste0("Finding correlation for ", rownames(act_CPM[i,]),
               " and mean_L1_CPM"))
  
  # create data frame with activator CPMs and mean L1 CPMs:
  df <- data.frame(as.numeric(act_CPM[i,]), as.numeric(L1_CPM_mean))
  colnames(df) <- c("activator_CPM", "mean_L1_CPM")
  df <- df[apply(df, 1, function(row) all(row !=0 )),]
  
  # calculate Spearman correlation:
  act_cor <- cor(df, method = "spearman")[2,1]
  
  if (i==1) {
    act_cors <- list(list(df, act_cor))
    names(act_cors)[i] <- rownames(act_CPM[i,])
  } else {
    act_cors[[i]] <- list(df, act_cor)
    names(act_cors)[i] <- rownames(act_CPM[i,])
  }
  
  # determine symbol corresponding to ID:
  m <- match(rownames(act_CPM[i,]), egENSEMBL$ensembl_id)
  gene_id <- egENSEMBL$gene_id[m]
  m <- match(gene_id, egSYMBOL$gene_id)
  sym <- egSYMBOL$symbol[m]
  
  # if correlation > 0.3, save as scatterplot:
  if ( act_cor > 0.3 ) {
    log_df <- log10(df)
    
    reg <- lm(mean_L1_CPM~activator_CPM, data = log_df)
    
    # reg <- lm(activator_CPM~mean_L1_CPM, data = log_df)
    xmax <- max(log_df$activator_CPM)
    ymax <- max(log_df$mean_L1_CPM)
    
    pdf(paste0(plotDir, "/", sym, 
               "_activator_correlation_with_mean_L1_CPM_exp.pdf"))
    par(mar = c(5,5,3,3))
    plot(log_df, ann = F)
    abline(reg, col="red")
    text(x=xmax-(0.1), y=ymax-(0.1), paste0("R=", round(act_cor, 3)))
    title(main=paste0("mean L1 vs ", sym), xlab = paste0(sym, 
                                                         " activator CPM"), ylab = paste0("mean L1 CPM"))
    dev.off()
    
    if ( !(exists("correlated_acts_0.3")) ) {
      correlated_acts_0.3 <- c(paste0("mean_L1_vs_", sym))
    } else {
      correlated_acts_0.3 <- c(correlated_acts_0.3, paste0("mean_L1_vs_", sym))
    }
  }
  
  # if correlation > 0.5, save as table:
  if ( act_cor > 0.5 ) {
    if ( !(exists("correlated_acts_0.5")) ) {
      correlated_acts_0.5 <- c(paste0("mean_L1_vs_", sym))
    } else {
      correlated_acts_0.5 <- c(correlated_acts_0.5, paste0("mean_L1_vs_", sym))
    }
  }
}
write.table(correlated_acts_0.3, paste0(tableDir, 
                                        "/L1_activators_correlated_with_mean_L1_CPMs_0.3.txt"), quote = F, sep = "\n")
write.table(correlated_acts_0.5, paste0(tableDir, 
                                        "/L1_activators_correlated_with_mean_L1_CPMs_0.5.txt"), quote = F, sep = "\n")


ctl_genes <- c("ENSG00000075624", "ENSG00000169919", "ENSG00000111640",
               "ENSG00000105173")
ctl_CPM <- CPM[as.character(ctl_genes),]
ctl_CPM <- na.omit(ctl_CPM)

for ( i in 1:nrow(ctl_CPM) ) {
  print(paste0("Finding correlation for ", rownames(ctl_CPM[i,]),
               " and mean L1 CPM"))
  
  # create data frame with ctl CPMs and mean L1 CPMs:
  df <- data.frame(as.numeric(ctl_CPM[i,]), as.numeric(L1_CPM_mean))
  colnames(df) <- c("ctl_CPM", "L1_CPM")
  df <- df[apply(df, 1, function(row) all(row !=0 )),]
  
  # calculate Spearman correlation:
  ctl_cor <- cor(df, method = "spearman")[2,1]
  
  if (i==1) {
    ctl_cors <- list(list(df, ctl_cor))
    names(ctl_cors)[i] <- rownames(ctl_CPM[i,])
  } else {
    ctl_cors[[i]] <- list(df, ctl_cor)
    names(ctl_cors)[i] <- rownames(ctl_CPM[i,])
  }
  
  # determine symbol corresponding to ID:
  m <- match(rownames(ctl_CPM[i,]), egENSEMBL$ensembl_id)
  gene_id <- egENSEMBL$gene_id[m]
  m <- match(gene_id, egSYMBOL$gene_id)
  sym <- egSYMBOL$symbol[m]
  
  # save as scatterplot:
  log_df <- log10(df)
  
  reg <- lm(L1_CPM~ctl_CPM, data = log_df)
  
  # reg <- lm(ctl_CPM~mean_L1_CPM, data = log_df)
  xmax <- max(log_df$ctl_CPM)
  ymax <- max(log_df$L1_CPM)
  
  # pdf(paste0(plotDir, "/", sym, "_activator_correlation_with_", 
  #            "mean_L1_CPM_exp.pdf"))
  # par(mar = c(5,5,3,3))
  plot(log_df, ann = F)
  # abline(reg, col="red")
  # text(x=xmax-(0.1), y=ymax-(0.1), paste0("r=", round(ctl_cor, 3)))
  # title(main=paste0(L1_ID, " vs ", sym), xlab = paste0(sym, " suppressor CPM"),
  #       ylab = paste0(L1_ID, " CPM"))
  # dev.off()
  
  # if correlation > 0.5, save as table:
  if ( ctl_cor > 0.5 ) {
    if ( !(exists("correlated_ctls_0.5")) ) {
      correlated_ctls_0.5 <- c(paste0(L1_ID, "_vs_", sym))
    } else {
      correlated_ctls_0.5 <- c(correlated_ctls_0.5, paste0(L1_ID, "_vs_", sym))
    }
  }
  
  # if correlation > 0.3, save as table:
  if ( ctl_cor > 0.3 ) {
    if ( !(exists("correlated_ctls_0.3")) ) {
      correlated_ctls_0.3 <- c(paste0(L1_ID, "_vs_", sym))
    } else {
      correlated_ctls_0.3 <- c(correlated_ctls_0.3, paste0(L1_ID, "_vs_", sym))
    }
  }  
}
write.table(correlated_ctls_0.3, paste0(tableDir, 
                                        "/L1_ctls_correlated_with_L1_CPMs_0.3.txt"), quote = F, sep = "\n")
write.table(correlated_ctls_0.5, paste0(tableDir, 
                                        "/L1_ctls_correlated_with_L1_CPMs_0.5.txt"), quote = F, sep = "\n")



