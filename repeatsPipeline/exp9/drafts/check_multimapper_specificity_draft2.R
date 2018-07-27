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

# create separate L1MD3 annotation:
L1MD3 <- reduce(reps[reps$type == "L1MD3"])

# remove L1MD3 from reps:
reps <- reps[reps$type != "L1MD3"]

### PLEASE REDO ABOVE AND SAVE FRESH ###
save.image(paste0(Robject_dir,
 "/check_multimapper_specificity_image1.RData"))


### 1. Isolate L1MD3 reads and check how many other 
# repeats each maps to ###

# identify reads from bam whose best alignments include L1MD3:
L1MD3_hits <- findOverlaps(bam_gr, L1MD3)
L1MD3_reads <- bam_gr[queryHits(L1MD3_hits)]

# identify all reads from bam_gr with same read names as in 
# L1MD3_reads:
alt_reads <- bam_gr[bam_gr$read_name %in% L1MD3_reads$read_name]

# remove L1MD3 regions and split into list by read name:
alt_reads  <- alt_reads[!(start(ranges(alt_reads)) %in% 
	start(ranges(L1MD3_reads)))]
alt_read_list <- split(alt_reads, alt_reads$read_name)

# view number of reads with multiple best alignments and number 
# of alternate best alignments for each read:
lens <- unlist(lapply(alt_read_list, length))
names(lens) <- NULL
lens

# count how many reads align to L1MD3 in total
L1MD3_split <- split(L1MD3_reads, L1MD3_reads$read_name)
length(L1MD3_split)

# identify what these reads that map to both L1MD3 and other 
# regions also map to:
other_aligns <- lapply(alt_read_list, function(x) {
	hits <- findOverlaps(x, reps)
	return(reps[subjectHits(hits)])
})

# check what other repeats can be mapped to by reads with best 
# alignments:
other_align_classes <- unique(
	unlist(lapply(other_aligns, function(x) {
		return(unique(x$type))
	}))
)
names(other_align_classes) <- NULL
other_align_classes

save.image(paste0(Robject_dir,
 "/check_multimapper_specificity_image2.RData"))


