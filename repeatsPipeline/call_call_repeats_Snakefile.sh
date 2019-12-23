#!/bin/bash

project_name="hgsoc_repeats"
exp_name="RNA-seq"

home_dir="/share/ScratchGeneral/jamtor"
project_dir="$home_dir/projects/$project_name/$exp_name"

 qsub -N repeats -b y -wd $project_dir -j y -R y -pe smp 44 \
 -V $project_dir/scripts/repeatsPipeline/call_repeats_Snakefile.sh
