#!/bin/bash

source activate p3.6.3env

project_name="hgsoc_repeats"
exp_name="RNA-seq"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$exp_name"

# call repeats snakemake:
snakemake -s $project_dir/repeats_Snakefile -k --cores 44
