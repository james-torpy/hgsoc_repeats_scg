
library(RColorBrewer)

htseq_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp9/htseq_EdgeR_primary_HGSOC_vs_FT_with_LINE1_silencers_extended/HGSOC_vs_FT_lrt.rds")

htseq_genes <- as.data.frame(topTags(htseq_lrt, n=Inf))
htseq_reps <- htseq_genes[grep("ENS",  rownames(htseq_genes), invert = T),]

all_reps <- list(htseq_reps)


salmon_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/salmonTE/SalmonTE_HGSOC_vs_FT/HGSOC_vs_FT_lrt.rds")

salmon_genes <- as.data.frame(topTags(salmon_lrt, n=Inf))
salmon_reps <- salmon_genes[grep("ENS",  rownames(salmon_genes), invert = T),]

all_reps[[2]] <- salmon_reps

all_reps[[3]] <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/TEtranscripts/custom3_DE.rds")

names(all_reps) <- c("htseq", "salmon", "TE")


rownames(all_reps[[3]]) <- gsub("\\:.*", "", rownames(all_reps[[3]]))

all_reps <- lapply(all_reps, function(x) {
  colnames(x) <- gsub("FDR|padj", "FDR_or_padj", colnames(x))
  return(subset(x, select = c(logFC, FDR_or_padj)))
})

plot_corr <- function(x, y) {
  df <- merge(x, y, by="row.names")
  rownames(df) <- df$Row.names
  df <- df[,-1]
  df$sig <- NA
  df$sig[df$FDR_or_padj.x > 0.1 & df$FDR_or_padj.y > 0.1] <- "non-significant"
  df$sig[df$FDR_or_padj.x < 0.1 & df$FDR_or_padj.y > 0.1] <- "x_significant"
  df$sig[df$FDR_or_padj.x > 0.1 & df$FDR_or_padj.y < 0.1] <- "y-significant"
  df$sig[df$FDR_or_padj.x < 0.1 & df$FDR_or_padj.y < 0.1] <- "both-significant"
  df$sig <- factor(df$sig, levels = c("non-significant", "x_significant",
                                      "y-significant", "both-significant"))
  
  logFC <- df[,grep("logFC", colnames(df))]
  corr <- cor(logFC, method="spearman")
  
  p <- ggplot(df, aes(logFC.x, logFC.y, colour = sig))
  p <- p + geom_point()
  p <- p + scale_colour_manual(values = c("grey", brewer.pal(3,"Paired")[3], 
                                          "gold1", "blueviolet"))
  p <- p + annotate("text", x=0, y=5, parse=TRUE, 
                    label=paste0("Spearman", 
                    as.character(round(corr[1,2], 3))))
  
  return(list(p, df))
}

par(mar=c(2,2,2,2))
htseq_vs_salmon <- plot_corr(all_reps[[1]], all_reps[[2]])


htseq_vs_TEtranscripts <- merge(x, all_reps[[3]], by="row.names")
rownames(htseq_vs_TEtranscripts) <- htseq_vs_TEtranscripts$Row.names
htseq_vs_TEtranscripts <- htseq_vs_TEtranscripts[,-1]
colnames(htseq_vs_TEtranscripts) <- c("htseq_logFC", "htseq_FDR", "TEtranscripts_logFC", 
                               "TEtranscripts_FDR")
h_vs_t_logFC <- htseq_vs_TEtranscripts[,grep("logFC", colnames(htseq_vs_TEtranscripts))]
h_vs_t_corr <- cor(h_vs_t_logFC, method="spearman")
plot(h_vs_t_logFC)

salmon_vs_TEtranscripts <- merge(all_reps[[2]], all_reps[[3]], by="row.names")
rownames(salmon_vs_TEtranscripts) <- salmon_vs_TEtranscripts$Row.names
salmon_vs_TEtranscripts <- salmon_vs_TEtranscripts[,-1]
colnames(salmon_vs_TEtranscripts) <- c("salmon_logFC", "salmon_FDR", 
                                       "TEtranscripts_logFC", "TEtranscripts_FDR")
s_vs_t_logFC <- salmon_vs_TEtranscripts[,grep("logFC", 
                                              colnames(salmon_vs_TEtranscripts))]
s_vs_t_corr <- cor(s_vs_t_logFC, method="spearman")
plot(s_vs_t_logFC)
