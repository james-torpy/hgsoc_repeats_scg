library(Rsamtools)
library(rtracklayer)
library(GenomicRanges)

#home_dir <- "/Users/jamestorpy/clusterHome/"
home_dir <- "/share/ScratchGeneral/jamtor/"
project_dir <- paste0(home_dir, "/projects/hgsoc_repeats/RNA-seq/")
results_dir <- paste0(project_dir, "results/")
ref_dir <- paste0(project_dir, "/refs/repeats/")
Robject_dir <- paste0(project_dir, "/Robjects/")

in_dir <- paste0(results_dir, "/star/GC/exp9/all_mmappers/prPT8/")


### 0. Load aligned bam file containing multimapped alignments ###

what <- c("qname", "flag", "rname", "strand", "pos", "qwidth")
param <- ScanBamParam(what = what)

bamFile <- paste0(in_dir, "/Aligned.out.bam")
bam <- scanBam(bamFile, param=param)[[1]]

# convert bam to GRanges object:
bam_gr <- GRanges(
	seqname=bam$qname,
	ranges=IRanges(start=bam$pos, width=bam$qwidth),
	strand=bam$strand,
	rname=bam$rname,
	flag=bam$flag
)
# remove odd sequences from gr:
bam_gr <- bam_gr[grep("M|K|G", bam_gr$seqname, invert=T)]

# load in repeats annotation:
reps <- import(paste0(project_dir, "/human-89.repeats.gtf"))
# remove odd sequences from reps:
reps <- reps[grep("M|K|G", reps$seqname, invert=T)]

save.image(paste0(Robject_dir, "/check_multimapper_specificity_image1.RData"))

# create separate L1MD3 annotation:
L1MD3 <- reduce(reps[reps$ID == "L1MD3"])


### 1. Isolate L1MD3 reads and check how many other repeats each maps to ###

# identify reads from bam whose best alignments include L1MD3:
L1MD3_hits <- findOverlaps(bam_gr, L1MD3)
L1MD3_reads <- bam_gr[queryHits(L1MD3_hits)]

# identify L1MD3 reads which also map to other repeats:
other_hits <- findOverlaps(L1MD3_reads, reps)
other_reads <- L1MD3_reads[queryHits(other_hits)]


