library(dplyr)

strv <- read.csv(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/paper_supps/supp5.csv", header = T, fill = T)
for (i in 2:5) {
  print(i)
  if (i==2) {
    strv_df <- data.frame(as.character(strv[,i][3:116]))
    colnames(strv_df) <- as.character(strv[,i][2:116][1])
  } else {
    strv_df <- cbind(strv_df, data.frame(as.character(strv[,i][3:116])))
    colnames(strv_df)[i-1] <- as.character(strv[,i][2:116][1])
  }
}

j=5
for (i in 19:27) {
  print(i)
  strv_df <- cbind(strv_df, data.frame(as.numeric(strv[,i][3:116])))
  colnames(strv_df)[j] <- as.character(strv[,i][2:116][1])
  j <<- j+1
}

strv <- strv_df

colnames(strv) <- gsub(" ", "_", colnames(strv))

strv <- dplyr::filter(strv, Sample_time_point == "primary")
strv$Case_ID <- as.character(strv$Case_ID)

sampleKey <- read.table(file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/sampleKey.txt", header = F, sep = " ", fill = T)[-c(1,2),3:4]
#sampleKey$V5 <- gsub("^.*_[0-9][0-9][0-9]_", "", sampleKey$V4)
#sampleKey <- sampleKey[grep("EXT", sampleKey$V5, invert = T),]
sampleKey$V4 <- gsub(
  "^.*nopd.", "", gsub(
    "_IC.*$|_EXTERNA.*$", "", sampleKey$V4
  )
)
colnames(sampleKey) <- c("id", "Case_ID")
sampleKey <- sampleKey[grep("rcAF|pAF|lrcT|msST", sampleKey$id, invert = T),]

strv <- merge(strv, sampleKey, by = "Case_ID")
strv$id <- gsub(
  "[0-9]$|[0-9][0-9]$", "", paste0(strv$Case_ID, "_", strv$id)
)
strv <- dplyr::select(strv, -c(Sample_time_point, Chemotherapy_response))


bp <- read.csv(file="/Users/jamestorpy/Documents/Garvan/phd/projects/hgsoc_repeats/paper_supps/supp6.csv", 
               header = T, fill = T)

for (i in 4:119) {
  print(i)
  if (i==4) {
    bp_df <- data.frame( as.character(bp[,i][2]), sum(as.numeric(as.character(bp[,i][3:25]))), 
                         sum(as.numeric(as.character(bp[,i][128:150]))) )
  } else {
    bp_df <- rbind(
      bp_df, data.frame( as.character(bp[,i][2]), sum(as.numeric(as.character(bp[,i][3:25]))), 
                  sum(as.numeric(as.character(bp[,i][128:150]))) )
    )
  }
}
colnames(bp_df) <- c("Case_ID", "Clustered_Breakpoints", "Breakpoint")

bp_df <- na.omit(bp_df)
bp_df$Tumour_sample_ID <- gsub(
  "_", "", gsub(
    "^.*_[0-9][0-9][0-9]_", "", bp_df$Case_ID
  )
)

bp_df$Case_ID <- gsub(
  "_ICGC.*$", "", bp_df$Case_ID
)


all_df <- merge(strv, bp_df, by = "Tumour_sample_ID")
rownames(all_df) <- all_df$id
all_df <- dplyr::select(all_df, -c(Tumour_sample_ID, Case_ID.x, id, Case_ID.y))
colnames(all_df)[1] <- "Total_SV"

write.table(all_df, file="/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/raw_files/fullsamples/bowtell_primary/SV_counts.txt",
            quote = F, sep = "\t")
