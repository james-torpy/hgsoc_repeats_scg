#!/bin/bash

echo -e
echo "### 5.htseq_all.bash ###"
echo -e
echo -e

source /home/jamtor/.bashrc

module load gi/boost/1.53.0
module load gi/zlib/1.2.8
module load phuluu/samtools/1.4

date

cores=$1
bamFile=$2
outFile=$3

# define gff file name:
gffFile="/home/jamtor/genomes/repeats/all_rep.final.gff"

# create outDir:
outDir=`echo ${outFile%/*}`/
mkdir -p $outDir

echo -e
echo This is the bamFile:
echo $bamFile
echo -e
echo This is the outFile:
echo $outFile
echo -e
echo This is the gffFile:
echo $gffFile
echo -e

#perform analysis on bam files using the reference genome
htseq_line="htseq-count -f bam -i ID -t exon $bamFile $gffFile >> $outFile"
echo This is the htseq_line:
echo $htseq_line
echo -e

htseq-count -f bam -i ID -t exon --stranded=no -a 0 -m intersection-strict $bamFile $gffFile >> $outFile

date