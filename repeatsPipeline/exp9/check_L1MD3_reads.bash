
# isolate column 1 of count assignment  info from all_reps.out.sam
# convert bam to sam
# add column from above to sam:
awk '{print $1}' ./htseq/exp9/prPT8/sam_output/prPT8.custom3.out.sam > ./htseq/exp9/prPT8/sam_output/prPT8.custom3.out2.sam
samtools view ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.bam > ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.sam
pr -mts' ' ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.sam ./htseq/exp9/prPT8/sam_output/prPT8.custom3.out2.sam > ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.annot.sam


# isolate reads aligned to repeat regions:
cat ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.annot.sam | grep -vE "no_feature|not_unique" > ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.repeats.sam


# isolate reads aligned to L1MD3 and save to L1MD3_reads:
samtools view -H ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.bam > ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.L1MD3.sam
cat ./star/GC/exp9/prPT8/Aligned.novosortedByName.out.repeats.sam | grep "L1MD3" | awk '{print $1}' >> ./star/GC/exp9/prPT8/L1MD3_reads.txt

# remove files no longer needed:
rm Aligned.novosortedByName.out.annot.sam Aligned.novosortedByName.out.sam Aligned.novosortedByName.out.L1MD3.sam Aligned.novosortedByName.out.repeats.sam