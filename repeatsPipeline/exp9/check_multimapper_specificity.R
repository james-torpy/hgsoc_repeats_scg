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

# load in repeats annotation:
reps <- import(paste0(ref_dir, "/human-89.repeats.gtf"))
# remove odd sequences from reps:
reps <- reps[grep("M|K|G", seqnames(reps), invert=T)]
values(reps) <- subset(values(reps), select=type)

# create separate L1MD3 annotation:
L1MD3 <- reduce(reps[reps$type == "L1MD3"])

# remove L1MD3 from reps:
reps <- reps[reps$type != "L1MD3"]

# remove odd sequences from gr:
bam_gr <- bam_gr[grep("M|K|G", seqnames(bam_gr), invert=T)]

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

# remove L1MD3 regions:
alt_reads  <- alt_reads[!(start(ranges(alt_reads)) %in% 
	start(ranges(L1MD3_reads)))]

# identify what these reads that map to both L1MD3 and other 
# regions also map to:
other_hits <- findOverlaps(alt_reads, reps)
other_aligns <- reps[queryHits(other_hits)]

save.image(paste0(Robject_dir,
 "/check_multimapper_specificity_image2.RData"))


