#!/bin/bash

tool="htseq"

project_dir="/share/ScratchGeneral/jamtor/projects/hgsoc_repeats/RNA-seq"
in_dir="$project_dir/results/$tool"
rhino_dir="10.0.11.16:/bridge/Zilog/Zilog-Cancer-TumourProgression/jamtor"
out_dir="$rhino_dir/projects/hgsoc_repeats/RNA-seq/results/$tool"
nci_address="jt3341@raijin.nci.org.au"
nci_dir="/g/data1a/ku3/jt3341/projects/hgsoc_repeats/RNA-seq/results/$tool"
mdss_dir="jt3341/projects/hgsoc_repeats/RNA-seq/results/$tool"

cd $in_dir

for f in $in_dir/AOCS*; do 
	id=$(basename $f)

	echo -e
	echo "tarring $in_dir/$id"
	tar -zcvf "$in_dir/$id.tar.gz" "$in_dir/$id"

	#echo -e
	#echo "transferring $in_dir/$id.tar.gz to $out_dir"
	#rsync -avPS "./$id.tar.gz" $out_dir

	echo -e
	echo "transferring $in_dir/$id.tar.gz to $nci_address:$nci_dir"
	rsync -avPS "$in_dir/$id.tar.gz" "$nci_address:$nci_dir"

	echo -e
	echo "transferring $nci_address:$nci_dir/$id.tar.gz to mdss $mdss_dir"
	ssh $nci_address "mdss put $nci_dir/$id.tar.gz $mdss_dir"

	echo -e
	echo "removing $nci_dir/$id.tar.gz"
	ssh $nci_address "rm $nci_dir/$id.tar.gz"

	echo -e
	echo "removing $in_dir/$id.tar.gz"
	rm "$in_dir/$id.tar.gz"
	
 done
