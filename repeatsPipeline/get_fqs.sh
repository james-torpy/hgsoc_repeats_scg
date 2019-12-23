#!/bin/bash

todl=( AOCS-135-4 AOCS-131-2 AOCS-130-2 AOCS-125-2 AOCS-124-2 AOCS-123-2 AOCS-122-2 AOCS-095-4 AOCS-092-4 AOCS-092-2 AOCS-091-4 AOCS-086-4 AOCS-086-2 AOCS-065-4 AOCS-064-4 AOCS-060-2 AOCS-034-4 AOCS-004-2 )

nci_address="jt3341@raijin.nci.org.au"
nci_dir="jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastqs/"
fq_dir1="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/raw_files/fastqs/"
fq_dir2="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/raw_files"

for f in ${todl[@]}; do 
	echo $f
	ssh $nci_address "mdss get $nci_dir/$f\_1.fastq.gz $fq_dir1"
	rsync -avPS $nci_address:$fq_dir1/$f\_1.fastq.gz $fq_dir2
	ssh $nci_address "rm $fq_dir1/$f\_1.fastq.gz"

	ssh $nci_address "mdss get $nci_dir/$f\_2.fastq.gz $fq_dir1"
	rsync -avPS $nci_address:$fq_dir1/$f\_2.fastq.gz $fq_dir2
	ssh $nci_address "rm $fq_dir1/$f\_2.fastq.gz"
done;