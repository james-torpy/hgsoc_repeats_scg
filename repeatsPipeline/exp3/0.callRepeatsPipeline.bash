# specify cores/job types needed:
countcores=20
subType="short"

# make directory hierachy:
projectname="hgsoc_repeats"

homeDir="/share/ScratchGeneral/jamtor/"
projectDir="$homeDir/projects/$projectname/RNA-seq/"
bamDir="$projectDir/results/star/GC/comparison/genome_aligns/"
RobjectDir="$projectDir/Robjects"
RDir="/home/jamtor/local/lib/r/R-3.2.2/bin/R"

echo -e
echo This is the bamDir:
echo $bamDir

# scripts/logs directories:
scriptsPath="$projectDir/scripts/repeatsPipeline/exp3"
logDir="$scriptsPath/logs/"
mkdir -p $logDir

echo -e
echo "The logDir is:"
echo $logDir

# fetch input bamfiles:
directories=( $(ls $bamDir | grep -v sorted | grep -v subset) )

echo -e
echo "The inDirs are: "
echo ${directories[@]}
echo -e

# submit count jobs for each sample:
for inDir in ${directories[@]}; do

	# create unique IDs for each bamfile:
	uID=`basename $inDir`

	annot="c"
	# 5a.countOverlaps1.R:
	qsub -q $subType.q -N fCca_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V "$RDir --vanilla --args $uID $annot \
	< $scriptsPath/5a.countOverlaps1.R"
	# 5b.countOverlaps_c2.R:
	qsub -q $subType.q -N fCcb_$uID -b y -wd $logDir -j y -R y -P \
	GenomeInformatics -pe smp $countcores -V "$RDir --vanilla --args $uID $annot \
	$annot < $scriptsPath/5b.countOverlaps2.R"
	
#	annot="custom3"
#	# 5a.countOverlaps_c1.R:
#	qsub -q $subType.q -N fCcus_$uID -b y -wd $logDir -j y -R y -P \
#	GenomeInformatics -pe smp $countcores -V "$RDir --vanilla --args $uID \
#	< $scriptsPath/5a.countOverlaps1.R"
	# 5b.countOverlaps_c2.R:
#	qsub -q $subType.q -N iisc5b_$uID.a -hold_jid sc3b_$uID -b y -wd $logDir -j y -R y -P \
#	GenomeInformatics -pe smp $countcores -V "$RDir --vanilla --args $uID \
#	$annot < $scriptsPath/repeatsPipeline/exp2a/5b.countOverlaps_c2.R"
done;






