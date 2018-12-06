
library(RColorBrewer)
library(edgeR)
library(ggplot2)

plotDir <- "/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/method_compare/"
system(paste0("mkdir -p ", plotDir))

htseq_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/exp9/htseq_EdgeR_primary_HGSOC_vs_FT_with_LINE1_silencers_extended/HGSOC_vs_FT_lrt.rds")

htseq_genes <- as.data.frame(topTags(htseq_lrt, n=Inf))
htseq_reps <- htseq_genes[grep("ENS",  rownames(htseq_genes), invert = T),]

all_reps <- list(htseq_reps)


salmon_lrt <- readRDS("/Users/jamestorpy/clusterHome/projects/hgsoc_repeats/RNA-seq/Robjects/salmonTE/SalmonTE_HGSOC_vs_FT/HGSOC_vs_FT_lrt.rds")

salmon_genes <- as.data.frame(topTags(salmon_lrt, n=Inf))
salmon_reps <- salmon_genes[grep("ENS",  rownames(salmon_genes), invert = T),]

all_reps[[2]] <- salmon_reps

TEreps <- read.table("/Users/jamestorpy/clusterHome//projects/hgsoc_repeats/RNA-seq/results/TEtranscripts/primary_HGSOC_vs_FT//primary_HGSOC_vs_FT_gene_TE_analysis.txt",
  header=T)
TEreps <- TEreps[grep("ENSG", TEreps$id, invert = T),]
rownames(TEreps) <- TEreps$id

TEreps <- subset(TEreps, select = c(log2FoldChange, padj))
colnames(TEreps) <- c("logFC", "padj")

all_reps[[3]] <- TEreps

# create vector of simplified repeat symbols:
simp <-  gsub(
  "\\(|\\)", "", gsub(
    "\\:.*$", "", rownames(all_reps[[3]])
  )
)

# find repeat symbols that are duplicated when simplifying and change them 
# back to original symbol names to make them repGenes rownames:
dupes <- simp[duplicated(simp)]
for ( d in dupes ) {
  simp[grep(d, simp)] <- rownames(all_reps[[3]])[grep(d, rownames(all_reps[[3]]))]
}
rownames(all_reps[[3]]) <- simp


names(all_reps) <- c("htseq", "salmon", "TE")

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
  
  # fit linear model to x and y data for regression line:
  reg_line <- logFC[grep("non", df$sig, invert=T),]
  fit1 <- lm(reg_line$logFC.y~reg_line$logFC.x)
  
  
  p <- ggplot(df, aes(logFC.x, logFC.y, colour = sig))
  p <- p + geom_point()
  p <- p + scale_colour_manual(values = c("grey", brewer.pal(3,"Paired")[3], 
                                          "gold1", "blueviolet"))
  p <- p + annotate("text", x=0, y=2.5, parse=TRUE, 
                    label=paste0("Spearman", 
                                 as.character(round(corr[1,2], 3))))
  p <- p + stat_smooth(method = "lm", col = "black", se = F, size=0.5)
  return(list(p, df))
}

par(mar=c(2,2,2,2))
htseq_vs_salmon <- plot_corr(all_reps[[1]], all_reps[[2]])
pdf(paste0(plotDir, "htseq_vs_salmon_repeat_DE.pdf"))
print(htseq_vs_salmon[[1]])
dev.off()

htseq_vs_TEtranscripts <- plot_corr(all_reps[[1]], all_reps[[3]])
pdf(paste0(plotDir, "htseq_vs_TEtranscripts_DE.pdf"))
print(htseq_vs_TEtranscripts[[1]])
dev.off()

salmon_vs_TEtranscripts <- plot_corr(all_reps[[1]], all_reps[[3]])
pdf(paste0(plotDir, "salmon_vs_TEtranscripts_DE.pdf"))
print(salmon_vs_TEtranscripts[[1]])
dev.off()
