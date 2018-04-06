#!/bin/bash

### TEtranscripts.bash ###

# This script takes a control and non-control group and performs a
# differential expression analysis between them for both ENSG and repeat
# genes


### 0. Define parameters and directories ###

cores=8

project="hgsoc_repeats"
exp_name="TEtranscripts"
sub_exp_name="primary_HGSOC_vs_FT"
#sub_exp_name="CCNE_HRD_PT_HGSOC_vs_FT"
#sub_exp_name="CCNE_PT_vs_HRD_PT_3_each"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project/RNA-seq/"
results_dir="$project_dir/results/"
bam_dir="$results_dir/star/GC/exp9/"
out_dir="$results_dir/$exp_name/$sub_exp_name"
annot_dir="$project_dir/refs/TEtranscripts/"
log_dir="/share/ClusterShare/thingamajigs/jamtor/projects/hgsoc_repeats/RNA-seq/logs/$exp_name/"
ref_dir="$project_dir/refs/TEtranscripts/"

mkdir -p $out_dir
mkdir -p $log_dir

echo -e
echo "The bam_dir is: $bam_dir"
echo -e
echo "The out_dir is: $out_dir"
echo -e
echo "The annot_dir is: $annot_dir"
echo -e
echo "The log_dir is: $log_dir"


### 1. Fetch input files ###

FT_files=( $(ls $bam_dir/**/Aligned.novosortedByName.out.bam | grep "FT") )


HGSOC_files=( $(ls $bam_dir/**/Aligned.novosortedByName.out.bam | grep -E \
"arPT|erPT|prPT|rfPT" | grep -vE "prPT28") )

# fetch filenames of CCNEamp and HRD driven prPT samples:
#HGSOC_IDs=( "arPT1" "erPT1" "prPT20" "prPT1" "prPT2" "rfPT2" )


# fetch filenames of CCNEamp and HRD driven PT samples:
#HGSOC_IDs=( $(cat "$ref_dir/CCNE_or_HRD_PT_samples_my_id_only.txt" | \
#	grep -vE "prPT28") )

#for i in $(seq 0 $((${#HGSOC_IDs[@]}-1)) ); do
#	echo $i
#
#	if [ $i -eq 0 ]; then
#		HGSOC_files=( "$bam_dir/${HGSOC_IDs[$i]}/Aligned.novosortedByName.out.bam" )
#	else
#		HGSOC_files+=( "$bam_dir/${HGSOC_IDs[$i]}/Aligned.novosortedByName.out.bam" )
#	fi;
#done;


#HGSOC_files=( ${HGSOC_files[@]:1:10} )

gencode_annot="$annot_dir/gencode_v24_hg38_annotation.gtf"
TE_annot="$annot_dir/hg38_rmsk_TE.gtf"

echo -e
echo "The FT_files are:"
echo ${FT_files[@]}
echo -e
echo "The HGSOC_files are:"
echo ${HGSOC_files[@]}
echo -e
echo "The gencode annotation is: $gencode_annot"
echo -e
echo "The TE annotation is: $TE_annot"


### 2. Build command and submit to cluster ###

TE_line="/home/jamtor/local/lib/TEToolkit-1.5.1/bin/TEtranscripts -t \
${HGSOC_files[@]} -c ${FT_files[@]} --GTF $gencode_annot --TE $TE_annot 
--project $out_dir"

echo -e
echo "The TE_line is:"
echo $TE_line


# activate python 2.7:
echo -e
echo "Activating python 2.7"
source activate p2.7env


# submit the job to the cluster:
qsub -N PT_TEtranscripts -b y -wd $log_dir -j y -R y -pe smp $cores -V $TE_line


### 1. Fetch CCNE vs HRD input files ###

#HRD_IDs=( "prPT1" "prPT2" "rfPT2" )
##HRD_IDs=( "prPT1" )
#
#for i in $(seq 0 $((${#HRD_IDs[@]}-1)) ); do
#	echo $i
#
#	if [ $i -eq 0 ]; then
#		HRD_files=( "$bam_dir/${HRD_IDs[$i]}/Aligned.novosortedByName.out.bam" )
#	else
#		HRD_files+=( "$bam_dir/${HRD_IDs[$i]}/Aligned.novosortedByName.out.bam" )
#	fi;
#done;
#
#CCNE_IDs=( "arPT1" "erPT1" "prPT20" )
##CCNE_IDs=( "arPT1" )
#
#for i in $(seq 0 $((${#CCNE_IDs[@]}-1)) ); do
#	echo $i
#
#	if [ $i -eq 0 ]; then
#		CCNE_files=( "$bam_dir/${CCNE_IDs[$i]}/Aligned.novosortedByName.out.bam" )
#	else
#		CCNE_files+=( "$bam_dir/${CCNE_IDs[$i]}/Aligned.novosortedByName.out.bam" )
#	fi;
#done;
#
#gencode_annot="$annot_dir/gencode_v24_hg38_annotation.gtf"
#TE_annot="$annot_dir/hg38_rmsk_TE.gtf"
#
#echo -e
#echo "The HRD_files are:"
#echo ${HRD_files[@]}
#echo -e
#echo "The CCNE_files are:"
#echo ${CCNE_files[@]}
#echo -e
#echo "The gencode annotation is: $gencode_annot"
#echo -e
#echo "The TE annotation is: $TE_annot"
#
#
#### 2. Build command and submit to cluster ###
#
#TE_line="/home/jamtor/local/lib/TEToolkit-1.5.1/bin/TEtranscripts -t \
#${CCNE_files[@]} -c ${HRD_files[@]} --GTF $gencode_annot --TE $TE_annot 
#--project $out_dir"
#
#echo -e
#echo "The TE_line is:"
#echo $TE_line
#
## activate python 2.7:
#echo -e
#echo "Activating python 2.7"
#source activate p2.7env
#
## submit the job to the cluster:
#qsub -N CCNEvsHRD_3each_TEtranscripts -b y -wd $log_dir -j y -R y -pe smp $cores -V $TE_line










