### 20.expression_vs_str_var_cats.R ###

# This script takes DE information for each primary HGSOC sample vs 
# pooled FT controls and plots it against number of SVs
# of different categories

### 0. Define variables/paths ###

# load packages needed:
library(plyr)
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
descrip <- paste0("expression_vs_structural_rearrangement_",
    "cats_indivHGSOC_grouped_unknown_only")
Type <- "custom3"
metric <- "HGSOC_vs_FT"
#repeatList <- c("(CATTC)n", "HSATII", "L1M3a", "L1M3c", "L1M4c", "L1M5", "L1M8", "L1MA5", "L1MC")
repeatPattern <- ""
FDRthresh <- 0.05
corr_thresh <- 2
remove_samples <- c("pAF", "rcAF")

divide_by_group <- TRUE
unknown_only <- TRUE

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

# define function to summarize data:
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


### 1. Load in DESeq results and create logFC data frame ###

if ( !file.exists(paste0(RobjectDir, "/repeat_logFC.rds")) & 
     !file.exists(paste0(RobjectDir, "/repeat_padj.rds")) &
     !file.exists(paste0(RobjectDir, "/orders.rds")) ) {
  
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
  
  saveRDS(FC, file=paste0(RobjectDir, "/repeat_FC.rds"))  
  
  
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
  
  saveRDS(orders, file=paste0(RobjectDir, "/orders.rds"))
  
} else {
  
  # load structural variation counts:
  str_var <- read.table(file=paste0(rawDir, "/fullsamples/bowtell_primary/SV_counts.txt"),
                        header = T, fill = T, sep = "\t")
  
  # load other necessary data frames:
  logFC <- readRDS(file=paste0(RobjectDir, "/repeat_logFC.rds"))
  padj <- readRDS(file=paste0(RobjectDir, "/repeat_padj.rds"))
  FC <- readRDS(file=paste0(RobjectDir, "/repeat_FC.rds"))
  orders <- readRDS(file=paste0(RobjectDir, "/orders.rds"))
  
}


### 3. Calculate CPMRs and test for correlations ###

Counts <- readRDS(file=paste0(RobjectDir, "/", Type, 
  "_counts.rds"))

# remove unwanted samples in remove_samples:
for ( r in remove_samples ) {
  print(r)
  Counts <- Counts[,grep(r, colnames(Counts), invert = T)]
}


# remove all HGSOC samples with HRD or CCNE1 gain if
# necessary:
if ( unknown_only ) {
  rem <- read.table(paste0(rawDir, 
    "/fullsamples/bowtell_primary", 
    "/CCNE_or_HRD_PT_samples_AOCS_id_only.txt"))
  nam <- gsub("\\_[a-zA-Z].*$", "", colnames(Counts))
  Counts <- Counts[,!(nam %in% rem$V1)]
}

# calculate total repeat count size:
rSizes <- apply(Counts, 2, sum)

# calculate CPMRs
CPMR <- as.data.frame(t(t(Counts)/rSizes)*1000000)

# only keep primary HGSOC samples:
CPMR <- CPMR[, grep("PT", colnames(CPMR))]

# only keep repeat information:
CPMR <- CPMR[grep("ENSG", rownames(CPMR), invert = T),]

colnames(CPMR) <- paste0(colnames(CPMR), "_vs_FT")

# clear variables:
rm(DE_cor)
rm(DE_cor_DE_rep)
rm(CPMR_cor)
rm(CPMR_cor_DE_rep)

# unlist orders:

for (m in 1:length(orders)) {
  print(paste0("m is ", m))
  # for each structural feature, order the logFC, FC and padj data frames by samples with
  # least to most of that feature:
  IDs <- paste0(orders[[m]], "_vs_FT")
  mtch <- na.omit(match(IDs, colnames(logFC)))
  logFC <- logFC[,mtch]
  padj <- padj[,mtch]
  FC <- FC[,mtch]
  
  for (j in 1:nrow(logFC)) {
    rp <- as.data.frame(t(logFC[j,]))
    rp$sample <- factor(rownames(rp), levels = rownames(rp))
    rp$padj <- as.numeric(t(padj[j,]))
    rp$pThresh <- as.logical((rp$padj < 0.05))
    rp$sig <- NA
    rp$sig[rp$pThresh] <- "< 0.05"
    rp$sig[!(rp$pThresh)] <- "non-significant"
      
    rp$FC <- as.numeric(t(FC[j,]))
    rp$str_no <- str_var[gsub(
      "_vs_FT", "", rownames(rp)), ][m]
      
    colnames(rp) <- c("log2FC", "sample", "padj", "pThresh", "sig", "FC", "str_no")
    if (j==1) {
      rpList <- list(rp)
    } else {
      rpList[[j]] <- rp
    }
  }
  names(rpList) <- rownames(logFC)[1:length(rpList)]
    
  n=1
  corList <- lapply(rpList, function(x) {
    df <- data.frame(x$str_no, x$FC)
    n <<- n+1
    return(cor(df, method = "spearman")[1,2])
  })
  names(corList) <- names(rpList)
  cors <- unlist(corList)
  
  rpInt <- rpList[cors > 0.1 | cors < -0.1]
  
  # if ( length(rpInt) > 0 ) {
  #   print("rpInt has length 0")
  #   break()
  # }

  
  ### 4. Plot logFC vs structural rearrangement for each repeat ###
  
  # fetch list of DE repeats for primary HGSOC vs FT:
  DEreps <- readRDS(file = paste0(RobjectDir, "/DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/custom3_DEsigReps.rds"))
  DEreps <- rownames(DEreps[[1]])
  
  for (k in 1:length(rpInt)) {
    if ( length(rpInt) > 0 ) {
      if ( !is.null(rpInt[[k]]) ) {
        rpName <- gsub("\\(|\\)", "", names(rpInt)[k])
        if ( !file.exists(paste0(plotDir, "/", rpName, "_DE_", names(orders)[m], "_0.4thresh.pdf")) ) {
          pdf(file = paste0(plotDir, "/", rpName, "_DE_", names(orders)[m], "_0.4thresh.pdf"), width=10, height=10)
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
        
        # save all repeat names with DE correlated with structural rearrangements to data frame:
        if ( !exists("DE_cor") ) {
          DE_cor <- data.frame(names(rpInt)[k], names(orders)[m], cors[names(rpInt)[k]])
        } else {
          DE_cor[,1] <- as.character(DE_cor[,1])
          DE_cor[,2] <- as.character(DE_cor[,2])
          DE_cor[,3] <- as.character(DE_cor[,3])
          DE_cor[(nrow(DE_cor)+1),] <- rep(NA, ncol(DE_cor))
          DE_cor[nrow(DE_cor), 1] <- names(rpInt)[k]
          DE_cor[(nrow(DE_cor)), 2] <- names(orders)[m]
          DE_cor[(nrow(DE_cor)), 3] <- cors[names(rpInt)[k]]
        }

        # check if repeat is differentially expressed in primary HGSOC vs FT:
        if ( rpName %in% DEreps ) {
          if ( !exists("CPMR_cor_DE_rep") ) {
            DE_cor_DE_rep <- data.frame(rpName, names(orders)[m])
          } else {
            DE_cor_DE_rep[,1] <- as.character(DE_cor_DE_rep[,1])
            DE_cor_DE_rep[,2] <- as.character(DE_cor_DE_rep[,2])
            DE_cor_DE_rep[(nrow(DE_cor_DE_rep)+1),] <- c(NA, NA)
            DE_cor_DE_rep[nrow(DE_cor_DE_rep), 1] <- rpName
            DE_cor_DE_rep[(nrow(DE_cor_DE_rep)), 2] <- names(orders)[m]
          }
        }
      }
    }
  }

  
  ### 5. Repeat for CPMR ###
  
  mtch <- na.omit(match(IDs, colnames(CPMR)))
  CPMR_ordered <- CPMR[,mtch]

  for (j in 1:nrow(CPMR_ordered)) {
    rp <- as.data.frame(t(CPMR_ordered[j,]))
    rp$sample <- factor(rownames(rp), levels = rownames(rp))
    
    rp$str_no <- str_var[gsub("_vs_FT", "", rownames(rp)), ][m]
    
    colnames(rp) <- c("CPMR", "sample", "str_no")
    if (j==1) {
      rpList <- list(rp)
    } else {
      rpList[[j]] <- rp
    }
  }
  names(rpList) <- rownames(CPMR_ordered)[1:length(rpList)]
  
  n=1
  corList <- lapply(rpList, function(x) {
    print(n)
    df <- data.frame(x$str_no, x$CPMR)
    n <<- n+1
    return(cor(df, method = "spearman")[1,2])
  })
  names(corList) <- names(rpList)
  cors <- unlist(corList)
  
  rpInt <- rpList[cors > 0.1 | cors < -0.1]
  rpInt <- rpInt[!is.na(names(rpInt))]
  print(names(rpInt))
  
  if ( divide_by_group ) {
    # group samples, 3 per group:
    for ( z in 1:3 ) {
      if (z==1) {
        no_vec <- c(rep(z, 6))
      } else {
        no_vec <- c(no_vec, rep(z, 6))
      }
    }
    no_vec <- no_vec[1:17]
    
    rpInt <- lapply(rpInt, function(x) {
      x$no <- no_vec

      sp <- split(x, x$no)
      x_se <- lapply(sp, function(y) {
        y <- subset(y, select=c(CPMR, sample, no))
        summarySE(y, measurevar="CPMR")
      })
      res <- do.call("rbind", x_se)
      colnames(res) <- c("group", "N", "mean_CPMR", "sd", "se", "ci")
      res$group <- 1:3
      return(res)
    })
    
    ### 4. Plot CPMR vs structural rearrangement for each repeat ###
    
    for (k in 1:length(rpInt)) {
      if ( length(rpInt) > 0 ) {
        if ( !is.null(rpInt[[k]]) ) {
          print(k)
          rpName <- gsub("\\(|\\)", "", names(rpInt)[k])
          if ( !file.exists(paste0(plotDir, "/", rpName, "_CPMR_", names(orders)[m], ".pdf"))) {
            pdf(file = paste0(plotDir, "/", rpName, "_CPMR_", names(orders)[m], ".pdf"), width=10, height=10)
            p <- ggplot(rpInt[[k]], aes(x = rpInt[[k]]$group, y = rpInt[[k]]$mean_CPMR))
            p <- p + geom_bar(stat="identity", position = "dodge")
            p <- p + geom_errorbar(aes(ymin=mean_CPMR-se, ymax=mean_CPMR+se),
                                   width=.2,                    # Width of the error bars
                                   position=position_dodge(.9))
            p <- p + ylab("mean CPMR vs FT")
            p <- p + xlab("primary HGSOC sample group (n=6)")
            p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
            p <- p + annotate("text", x=5, y=300, label=paste0("Spearman = ", as.character(round(cors[rpName], 4))))
            p <- p + ggtitle(paste0(rpName, " vs ", names(orders)[m]))
            print(p)
            dev.off()
          }
        }
      }
    }  
  }
  
  
  
  ### 4. Plot CPMR vs structural rearrangement for each repeat ###

  for (k in 1:length(rpInt)) {
    if ( length(rpInt) > 0 ) {
      if ( !is.null(rpInt[[k]]) ) {
        print(k)
        rpName <- gsub("\\(|\\)", "", names(rpInt)[k])
        if ( !file.exists(paste0(plotDir, "/", rpName, "_CPMR_", names(orders)[m], ".pdf"))) {
          pdf(file = paste0(plotDir, "/", rpName, "_CPMR_", names(orders)[m], ".pdf"), width=10, height=10)
          p <- ggplot(rpInt[[k]], aes(x = rpInt[[k]]$sample, y = rpInt[[k]]$CPMR))
          p <- p + geom_bar(stat="identity", position = "dodge")
          p <- p + ylab("CPMR vs FT")
          p <- p + xlab("sample")
          p <- p + theme(axis.text.x = element_text(angle = 90, hjust = 1))
          p <- p + annotate("text", x=45, y=7, label=paste0("Spearman = ", as.character(cors[rpName])))
          p <- p + ggtitle(paste0(rpName, " vs ", names(orders)[m]))
          print(p)
          dev.off()
        }
      
        # save all repeat names with CPMR correlated with structural rearrangements to data frame:
        if ( !exists("CPMR_cor") ) {
          CPMR_cor <- data.frame(names(rpInt)[k], names(orders)[m], cors[names(rpInt)[k]])
        } else {
          CPMR_cor[,1] <- as.character(CPMR_cor[,1])
          CPMR_cor[,2] <- as.character(CPMR_cor[,2])
          CPMR_cor[,3] <- as.character(CPMR_cor[,3])
          CPMR_cor[(nrow(CPMR_cor)+1),] <- rep(NA, ncol(CPMR_cor))
          CPMR_cor[nrow(CPMR_cor), 1] <- names(rpInt)[k]
          CPMR_cor[(nrow(CPMR_cor)), 2] <- names(orders)[m]
          CPMR_cor[(nrow(CPMR_cor)), 3] <- cors[names(rpInt)[k]]
        }
        
        # check if repeat is differentially expressed in primary HGSOC vs FT:
        if ( rpName %in% DEreps ) {
          if ( !exists("CPMR_cor_DE_rep") ) {
            CPMR_cor_DE_rep <- data.frame(rpName, names(orders)[m])
          } else {
            CPMR_cor_DE_rep[,1] <- as.character(CPMR_cor_DE_rep[,1])
            CPMR_cor_DE_rep[,2] <- as.character(CPMR_cor_DE_rep[,2])
            CPMR_cor_DE_rep[(nrow(CPMR_cor_DE_rep)+1),] <- c(NA, NA)
            CPMR_cor_DE_rep[nrow(CPMR_cor_DE_rep), 1] <- rpName
            CPMR_cor_DE_rep[(nrow(CPMR_cor_DE_rep)), 2] <- names(orders)[m]
          }
        }
      }
    }
  }
}

print("Repeats with correlation between structural rearrangements and DE values and showing DE for HGSOC vs FT are:")
DE_cor_DE_rep

print("Repeats with correlation between structural rearrangements and CPMR values and showing DE for HGSOC vs FT are:")
CPMR_cor_DE_rep

# save DE_cor and CPMR_cor as tab-delimited files:
colnames(DE_cor) <- c("repeat_id", "rearrangement", "correlation with expression")
if ( !file.exists(paste0(plotDir, "/repeats_DE_cor_with_str_var.txt")) ) {
  write.table(DE_cor, file=paste0(plotDir, "/repeats_DE_cor_with_str_var.txt"), quote = F, sep = "\t", 
            row.names = F, col.names = F)
} else {
  print(paste0(plotDir, "/repeats_DE_cor_with_str_var.txt", " already exists"))
}

colnames(CPMR_cor) <- c("repeat_id", "rearrangement", "correlation with expression")
if ( !file.exists(paste0(plotDir, "/repeats_CPMR_cor_with_str_var.txt")) ) {
  write.table(CPMR_cor, file=paste0(plotDir, "/repeats_CPMR_cor_with_str_var.txt"), quote = F, sep = "\t", 
            row.names = F, col.names = F)
} else {
  print(paste0(plotDir, "/repeats_CPMR_cor_with_str_var.txt", " already exists"))
}

# fetch FT vs HGSOC DE data for each correlated repeat:
DE <- readRDS(paste0(projectDir, "/RNA-seq/Robjects/exp9/DEverify/htseq_EdgeR_primary_HGSOC_vs_FT/", 
  "custom3_DEreps.rds"))[[1]]

CPMR_cor$logFC <- DE[CPMR_cor$repeat_id,]$logFC
CPMR_cor$FDR <- DE[CPMR_cor$repeat_id,]$FDR

saveRDS(CPMR_cor, paste0(RobjectDir, "/CPMR_cor.rds"))
