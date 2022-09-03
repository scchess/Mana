# Test bash shell for complete binf pipeline

# Set the path where this bash shell will be
shellpath=/home/uqsidris/'v2Pipeline_plasmid_200722.sh'

# Set the path to look in for passed fastq files
passpath=/home/uqsidris/DATA_RUNS/Project_017_LSK109_plasmids/mRNA14_plasmid_mRNA17_plasmid/20220525_0136_X2_AJS173_301c748a/fastq_pass

# Set the path to look in for failed fastq files
failpath=/home/uqsidris/DATA_RUNS/Project_017_LSK109_plasmids/mRNA14_plasmid_mRNA17_plasmid/20220525_0136_X2_AJS173_301c748a/fastq_fail

# Set the path to save outputs
outputpath=/home/uqsidris/workdir/plasmid_flongle

# Set the path to look in for the reference fasta files
referencepath=/home/uqsidris/workdir/referencefiles

# Set path to seq run summary file
runsummaryfilepath=/home/uqsidris/DATA_RUNS/Project_017_LSK109_plasmids/mRNA14_plasmid_mRNA17_plasmid/20220525_0136_X2_AJS173_301c748a/sequencing_summary_AJS173_1fa1ba1d

# Set array for the Barcodes used match format to the name of the folders that the data is in
barcode_arr=(barcode10 barcode11)

# Set array for the sample names that match in the order of the barcodes in the array above
sample_arr=(mrna14_plasmid mrna17_plasmid)

# Set array for the reference files in fasta format that match in the order of the samples in the array above
reference_arr=(2022_14_482162_F4 2022_17_483289_F8)

# Create folders with sample names in the output directory
for sample in ${sample_arr[@]}; do
    mkdir -p $outputpath/$sample
done

# make logfile
LOG_FILE=$outputpath/audit_trail.log
exec > >(while read -r line; do printf '%s %s\n' "$(date --rfc-3339=seconds)" "$line" | tee -a $LOG_FILE; done)
exec 2> >(while read -r line; do printf '%s %s\n' "$(date --rfc-3339=seconds)" "$line" | tee -a $LOG_FILE; done >&2)

# copy the shell used as a record
cp $shellpath $outputpath

# nanoplot for the whole run and each sample
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    NanoPlot --summary $runsummaryfilepath'.txt' --loglength -o $outputpath/$samplefile/$samplefile'run_nanoplot'
done

# Join passed fastq files from each barcode then rename and save in appropriate sample folder
for i in ${!barcode_arr[@]}; do
    barcodefile=${barcode_arr[i]}
    samplefile=${sample_arr[i]}
    zcat $passpath/$barcodefile/*fastq.gz > $outputpath/$samplefile/$samplefile'_beforetrimming_allpassedreads.fastq.gz'
done

# Join failed fastq files from each barcode then rename and save in appropriate sample folder
for i in ${!barcode_arr[@]}; do
    barcodefile=${barcode_arr[i]}
    samplefile=${sample_arr[i]}
    zcat $failpath/$barcodefile/*fastq.gz > $outputpath/$samplefile/$samplefile'_allfailedreads.fastq.gz'
done

porechop --check_reads 5000000 -i $outputpath/$samplefile/$samplefile'_beforetrimming_allpassedreads.fastq.gz' -o $outputpath/$samplefile/$samplefile'_allpassedreads.fastq.gz'

# Running minimap2 align fastq to reference fasta and make sam file in correct directory for all passed reads
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    minimap2 -ax map-ont $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads.fastq.gz' > $outputpath/$samplefile/$samplefile'_allpassedreads.sam'
done

# Running minimap2 align fastq to reference fasta and make sam file in correct directory for all failed reads
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    minimap2 -ax map-ont $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allfailedreads.fastq.gz' > $outputpath/$samplefile/$samplefile'_allfailedreads.sam'
done

-----------------------------------------------------------
                    MANA STARTS HERE
-----------------------------------------------------------

# Run samtools for various statistics and also create a sorted and indexed bam file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # create samtools stats file for the sam file for all passed reads
    samtools stats $outputpath/$samplefile/$samplefile'_allpassedreads.sam' > $outputpath/$samplefile/$samplefile'_allpassedreads_stats.txt'
    samtools stats $outputpath/$samplefile/$samplefile'_allfailedreads.sam' > $outputpath/$samplefile/$samplefile'_allfailedreads_stats.txt'
    # create a bam file from the sam file for all passed reads
    samtools view -S -b $outputpath/$samplefile/$samplefile'_allpassedreads.sam' > $outputpath/$samplefile/$samplefile'_allpassedreads.bam'
    samtools view -S -b $outputpath/$samplefile/$samplefile'_allfailedreads.sam' > $outputpath/$samplefile/$samplefile'_allfailedreads.bam'
    # sort the bam file for all passed reads
    samtools sort $outputpath/$samplefile/$samplefile'_allpassedreads.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam'
    samtools sort $outputpath/$samplefile/$samplefile'_allfailedreads.bam' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam'
    # index the bam file for all passed reads
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam'
    # create samtools flagstat file for all passed reads
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat.txt'
done

# Run pysamstats to get various errors per base then use awk to remove obsolete columns and create new file then do some calculations and output a new file
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    # make pysamstats txt file
    pysamstats --max-depth=300000000 --fasta $referencepath/$referencefile'.fasta' --type variation $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_pysam.txt'
    # all together cut and percent each error in new column with headers
    cut -f 5,7,9,11,13,15,17,19,21,23 --complement $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR==1{$13 = $13 OFS "Match per base" OFS "Mismatch per base" OFS "Deletion per base" OFS "Insertion per base"} NR>1{$13 = $13 OFS ($5/$4) OFS ($6/$4) OFS ($7/$4) OFS ($8/$4)}1' CONVFMT='%.3f' > $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt'
    # making a table with min max average and headers.
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$14}; if($14>max) {max=$14}; if($14<min) {min=$14}; total+=$14; count+=1} {CONVFMT="%.3f"} END {print "Error Performance per Nucleotide" "\n" "Frequency %" OFS "avg" OFS "max" OFS "min" "\n" "Match" OFS total/count OFS max OFS min}' > $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt'
    # Next Line next error insertion
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$17}; if($17>max) {max=$17}; if($17<min) {min=$17}; total+=$17; count+=1} {CONVFMT="%.3f"} END {print "Insertion" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt'
    # Next Line next error mismatch
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$15}; if($15>max) {max=$15}; if($15<min) {min=$15}; total+=$15; count+=1} {CONVFMT="%.3f"} END {print "Mismatch" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt'
    # Next Line next error deletion
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$16}; if($16>max) {max=$16}; if($16<min) {min=$16}; total+=$16; count+=1} {CONVFMT="%.3f"} END {print "Deletion" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt'
    # Next Line next total error
    cat $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR>4{if(min==""){min=$4}; if(max==""){max=$3}; if ($3>max) {max=$3}; if($4<min) {min=$4}; total+=$2; count+=1} {CONVFMT="%.3f"} END {print "Total Error" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt'
    # Notes
    cat $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt' | awk 'BEGIN{FS=OFS="\t"} END {print "Notes" "\n" "Error = mismatch + deletion" "\n" "Mismatch does not include deletion"}' >> $outputpath/$samplefile/$samplefile'_passed_pysamreporttable.txt'
done

# Run awk etc to get various numbers from flagstat, stats and pysamtable do some calculations and output a new file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # making a table with title Run performance
    echo "Run Performance" > $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    # convert flagstat file to tab delimited file
    cat $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat_tabdelimited.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat_tabdelimited.txt'
    # open tab delimited flagstat file and print total reads number (1st number in the file) for passed to the report
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$1} END {print "Passed Reads >Q9" OFS passedread}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{failedread=$1} END {print "Failed Reads <Q9" OFS failedread}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    # sum the above to get total reads from the run
    cat $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} NR==3{failedread=$2} END {print "Total Reads" OFS passedread+failedread}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    # find flagstat mapped number 7th number in file
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==8{mappedread=$1} END {print "Mapped Reads" OFS mappedread}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    # calc unmapped reads 
    cat $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} NR==5{mappedread=$2} END {print "Unmapped Reads" OFS passedread-mappedread}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    # convert samtools stats file to tab delimited file
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_stats.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_passed_stats_tabdelimited.txt'
    # find read length average from samtools stats
    cat $outputpath/$samplefile/$samplefile'_passed_stats_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==33{readlengthavg=$4} END {print "Read Length Average" OFS readlengthavg}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    # find coverage -average -max -min from the pysamstats report
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Average Coverage" OFS total/count}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Min Coverage" OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_passed_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Max Coverage" OFS max}' >> $outputpath/$samplefile/$samplefile'_passed_allreads_runperformance_reporttable.txt'
done

# Create file listing the lengths of each primary mapped read in the sorted bam file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools view -F 2048 $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' | awk '{print length($10)}'| sort -n | uniq -c |awk ' { t = $1; $1 = $2; $2 = t; print; } '|tr ' ' '\t'|sed '1d' > $outputpath/$samplefile/$samplefile'_passed_sorted_primarymappedreadlengths.txt'
    samtools view -F 2048 $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam' | awk '{print length($10)}'| sort -n | uniq -c |awk ' { t = $1; $1 = $2; $2 = t; print; } '|tr ' ' '\t'|sed '1d' > $outputpath/$samplefile/$samplefile'_failed_sorted_primarymappedreadlengths.txt'
done

# Make a bam of unmapped reads
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools view -S -b -f 4 $outputpath/$samplefile/$samplefile'_allpassedreads.sam' > $outputpath/$samplefile/$samplefile'_unmapped_passedreads.bam'
    samtools index $outputpath/$samplefile/$samplefile'_unmapped_passedreads.bam'
done

# Add variant calling and consensus generation
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' | bcftools call -c -M --ploidy 1 -Oz -o $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz'
    bcftools index $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz'
    bcftools norm -f $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz' -Ob -o $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_BCF.bcf'
    bcftools consensus -a - -f $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_BCF_consensus.fa'
done

