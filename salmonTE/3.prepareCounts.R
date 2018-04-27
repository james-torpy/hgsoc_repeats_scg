### 3.prepareCounts.R ###

# breaks gcCounts up into dataframe for each class and
# saves each as RData objects #

# Run on cluster with:
#briR
#qsub -N salPrep -b y -wd \
#/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/exp9/DE \
#-j y -R y -pe smp 4 -V "Rscript /share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/exp9/salmonTE/3.prepareCounts.R"


# define starting variables:
project <- "hgsoc_repeats"
expName <- "salmonTE"

# define directories:
#homeDir <- "/share/ScratchGeneral/jamtor/"
homeDir <- "/Users/jamestorpy/clusterHome"
projectDir <- paste0(homeDir, "/projects/", 
  project, "/")
resultsDir <- paste0(projectDir, "/RNA-seq/results/")
RobjectDir <- paste0(projectDir, "/RNA-seq/Robjects/", 
  expName, "/")
inDir <- paste0(resultsDir, "/", expName)
plotDir <- paste0(inDir, "/plots/")
rawDir <- paste0(projectDir, 
  "/RNA-seq/raw_files/fullsamples/bowtell_primary/")

system(paste0("mkdir -p ", RobjectDir))
system(paste0("mkdir -p ", plotDir))


### 1. Check for and discard bad inputs ###

# fetch counts file ids:
compFiles <- list.files(inDir, pattern = "[A-Z][A-Z}[0-9].*")

# remove non-primary samples:
compFiles <- grep("PT|FT", compFiles, value=T)

# select subset of files:
#compFiles <- grep(
#  "prPT[0-9]$", grep(
#    "FT3", compFiles, invert=T, value=T
#  ), invert=T, value=T
#)


### 2. Convert Salmon GC counts ###

if ( exists("CountsGC") ) {
  upTo <- ncol(CountsGC)-1
} else {
  upTo <- 1
}

if ( !file.exists(paste0(RobjectDir, "gcCounts_temp.rds")) ) {
  # load in gc_countFiles into a df:
  for (i in upTo:length(compFiles)) {
    if (i==1) {
      print(compFiles[i])
      CountsGC <- data.frame(read.table(file=paste0(inDir, "/", compFiles[i], "/hg38_transcriptome/quant.sf"), header=T))[,c(1,5)]
      colnames(CountsGC) <- c("gene_id", compFiles[i])
      print(compFiles[i])
    } else {
      print(compFiles[i])
      CountsGC <- cbind(CountsGC, read.table(file=paste0(inDir, "/", compFiles[i], "/hg38_transcriptome/quant.sf"), header=T)[,5])
      colnames(CountsGC)[i+1] <- compFiles[i]
      print(compFiles[i])
    }
  }
  
  CountsGC$gene_id <- gsub(
      "\\..*$", "", gsub(
        "\\|.*$", "", gsub(
        "^.*ENSG", "ENSG", CountsGC$gene_id
      )
    )
  )
  
  saveRDS(CountsGC, paste0(RobjectDir, "gcCounts_temp.rds"))
  
} else {
  
  CountsGC <- readRDS(paste0(RobjectDir, "gcCounts_temp.rds"))
}

Counts <- CountsGC
# aggregate multiple types of same gene:
Counts <- aggregate(.~gene_id, Counts, mean)
# remove duplicate rows:
Counts <- unique(Counts)

# replace uIDs with patient ids:
# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# renames duplicate sample names:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(Counts))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(Counts, select=-gene_id))))
newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
Key$V4 <- gsub("\\_[a-zA-Z].*$", "", newNames)

# change colnames of Counts:
m <- match(colnames(Counts), Key$V3)
colnames(Counts) <- gsub("NA_", "", 
                         gsub("[0-9]$", "",
                              gsub(
                                "[0-9][0-9]$", "", paste0(Key$V4[m], "_", colnames(Counts))
                              )
                         )
)


### 3. Convert SalmonTE counts ###

if ( exists("repCounts") ) {
  upTo <- ncol(repCounts)-1
} else {
  upTo <- 1
}

if ( !file.exists(paste0(RobjectDir, "SalmonTECounts_temp.rds")) ) {
  # load in gc_countFiles into a df:
  for (j in 1:length(compFiles)) {
    if (j==1) {
      print(compFiles[j])
      repCounts <- data.frame(read.csv(file=paste0(inDir, "/", compFiles[j], "/custom3/EXPR.csv"), header=T))
      colnames(repCounts) <- c("gene_id", compFiles[i])
    } else {
      print(compFiles[j])
      if ( file.exists(paste0(inDir, "/", compFiles[j], "/custom3/EXPR.csv")) ) {
        repCounts <- cbind(repCounts, read.csv(file=paste0(inDir, "/", compFiles[j], "/custom3/EXPR.csv"), header=T)[,2])
      } else {
        repCounts <- cbind(repCounts, rep(NA, nrow(repCounts)))
      }
      colnames(CountsGC)[i+1] <- compFiles[i]
    }
  }
  # name columns:
  colnames(repCounts) <- c("gene_id", compFiles)
  
  # omit NAs:
  repCounts <- as.data.frame(repCounts[,colSums(is.na(repCounts)) < nrow(repCounts)])
  
  # coerce data frame values to integers:
  repCounts <- as.data.frame(repCounts[,colSums(is.na(repCounts)) < nrow(repCounts)])
  nam <- repCounts$gene_id
  repCounts <- subset(repCounts, select=-gene_id)
  repCounts <- as.data.frame(sapply(repCounts, as.integer))
  repCounts$gene_id <- nam
  repCounts <- repCounts[,c( ncol(repCounts), seq(1, (ncol(repCounts)-1)) )]
  
  saveRDS(repCounts, paste0(RobjectDir, "SalmonTECounts_temp.rds"))
  
} else {
  
  repCounts <- readRDS(paste0(RobjectDir, "SalmonTECounts_temp.rds"))
}

# replace uIDs with patient ids:
# load in patient ids:
Key <- read.table(file=paste0(rawDir, "/sampleKey.txt"), sep=" ", fill=T)[,3:4]
Key <- Key[3:nrow(Key),]
Key$V4 <- gsub("_[A-Z].*$", "", gsub("^.*AOCS_", "AOCS_", Key$V4))
Key <- Key[with(Key, order(V3)), ]

# renames duplicate sample names:
Key <- dplyr::filter(Key, Key$V3 %in% colnames(repCounts))
newNames <- paste0(Key$V4, "_", gsub("[0-9]", "", colnames(subset(repCounts, select=-gene_id))))
newNames[duplicated(newNames)] <- gsub("AOCS", "AOCS2", newNames[duplicated(newNames)])
newNames[duplicated(newNames)] <- gsub("AOCS2", "AOCS3", newNames[duplicated(newNames)])
Key$V4 <- gsub("\\_[a-zA-Z].*$", "", newNames)

# change colnames of repCounts:
m <- match(colnames(repCounts), Key$V3)
colnames(repCounts) <- gsub("NA_", "", 
                         gsub("[0-9]$", "",
                              gsub(
                                "[0-9][0-9]$", "", paste0(Key$V4[m], "_", colnames(repCounts))
                              )
                         )
)

# remove some columns of Counts so it matches columns of repCounts:
Counts <- Counts[,colnames(Counts) %in% colnames(repCounts)]

# save Counts:
if (!file.exists(paste0(RobjectDir, "/gc_counts.Salmon.rds"))) {
  saveRDS(Counts, file=paste0(RobjectDir, "/gc_counts.Salmon.rds"))
}

# save repCounts:
if (!file.exists(paste0(RobjectDir, "/custom3_counts.SalmonTE.rds"))) {
  saveRDS(repCounts, file=paste0(RobjectDir, "/custom3_counts.SalmonTE.rds"))
}
