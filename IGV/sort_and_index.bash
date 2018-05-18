#!/bin/bash
#/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq/scripts/IGV/sort_and_index.bash

sample="FT7"

home_dir="/share/ScratchGeneral/jamtor/"
results_dir="$home_dir/projects/hgsoc_repeats/RNA-seq/results/"
in_dir="$results_dir/star/GC/exp9/$sample/"
out_dir="$results_dir/IGV/$sample/"
temp_dir="$out_dir/temp/"

mkdir -p $out_dir
mkdir -p $temp_dir

sort_line="/home/jamtor//local/bin/novocraft/novosort \
-t $temp_dir -c 6 -m 22G \
$in_dir/Aligned.novosortedByName.out.bam > \
$out_dir/Aligned.novosortedByCoordinate.out.bam; index \
$out_dir/Aligned.novosortedByCoordinate.out.bam"

echo $sort_line

qsub -N sort_$sample -b y -wd $out_dir -j y -R y -pe smp 6 -V $sort_line