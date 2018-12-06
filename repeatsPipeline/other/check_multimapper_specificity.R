library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

#home_dir <- "/Users/jamestorpy/clusterHome/"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, 
	"/projects/hgsoc_repeats/RNA-seq/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "/refs/repeats/")
Robject_dir <- paste0(project_dir, "/Robjects/")
in_dir <- paste0(results_dir, 
	"/star/GC/exp9/all_mmappers/prPT8/")
in_dir2 <- paste0(results_dir, "star/GC/exp9/prPT8")


### 0. Load aligned bam file containing multimapped alignments ###

what <- c("qname", "flag", "rname", "strand", "pos", "qwidth")
param <- ScanBamParam(what = what)

bamFile <- paste0(in_dir, "/Aligned.out.bam")
bam <- scanBam(bamFile, param=param)[[1]]

# convert bam to GRanges object:
bam_gr <- GRanges(
	seqnames=bam$rname,
	ranges=IRanges(start=bam$pos, width=bam$qwidth),
	strand=bam$strand,
	read_name=bam$qname,
	flag=bam$flag
)

# remove all entries marked as non-primary (flag = 256 or 272)
# as STAR was directed to mark all best score alignments as 
# primary (0 or 16):
bam_gr <- bam_gr[!(bam_gr$flag %in% c(256, 272))]

# remove odd sequences from gr:
bam_gr <- bam_gr[grep("M|K|G", seqnames(bam_gr), invert=T)]

# load in repeats annotation:
reps <- import(paste0(ref_dir, "/human-89.repeats.gtf"))
# remove odd sequences from reps:
reps <- reps[grep("M|K|G", seqnames(reps), invert=T)]
values(reps) <- subset(values(reps), select=type)
reps <- reps[!(reps$type %in% c("dust", "trf"))]

save.image(paste0(Robject_dir,
 "/check_multimapper_specificity_image1.RData"))


### 1. Load L1MD3 counted reads and check how many other 
# repeats each maps to ###

L1MD3_reads <- as.character(unique(read.table(paste0(in_dir2, 
	"/L1MD3_reads.txt"))[,1]))

# identify reads from bam whose best alignments include L1MD3:
alt_reads <- bam_gr[bam_gr$read_name %in% L1MD3_reads]

# split into list by read name:
alt_read_list <- split(alt_reads, alt_reads$read_name)

# identify what these reads that map to both L1MD3 and other 
# regions also map to:
other_aligns <- lapply(alt_read_list, function(x) {
	hits <- findOverlaps(x, reps)
	res <- reps[subjectHits(hits)]
	res$read_name <- x$read_name
	if (length(res) < 1) {
		res <- NULL
	}

	return(res)
})
other_aligns <- other_aligns[!sapply(other_aligns, is.null)]

minus_L1MD3 <- lapply(other_aligns, function(x) {
	if ( length(unique(x$type)) <= 1) {
		x <- NULL
	}
	return(x)
})
minus_L1MD3 <- minus_L1MD3[!sapply(minus_L1MD3, is.null)]

# check what other repeats can be mapped to by reads with best 
# alignments:
other_align_classes <- unique(
	unlist(lapply(other_aligns, function(x) {
		return(unique(x$type))
	}))
)
names(other_align_classes) <- NULL
other_align_classes

genes <- c("L1MC4a", "L1MD2", "L1MC3", "L1M4", "(GA)n", "L1MB5")

for ( i in 1:length(minus_L1MD3) ) {
	for ( j in 1:length(genes) ) {
		res <- grep(genes[j], minus_L1MD3[[i]]$type)
		if ( length(res) > 0 ) {
			print(unique(minus_L1MD3[[i]]$read_name))
			print(genes[j])
			print(res)
		}
	}
}


alt_read_list[grep("D81P8DQ1:153:C2704ACXX:3:2304:12859:29397", 
	names(alt_read_list))]

save.image(paste0(Robject_dir,
 "/check_multimapper_specificity_image2.RData"))