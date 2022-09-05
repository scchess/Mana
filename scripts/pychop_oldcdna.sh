# Test bash shell for complete binf pipeline

### Set the path where this bash shell will be
shellpath=/home/uqsidris/'pychop_oldcdna.sh'
### Set the path to look in for passed fastq files
passpath=/home/uqsidris/network/SEQLAB-Q4547/cDNA-PCR_mRNA_vaccine_QC_010/mRNA_vaccines/20220422_0640_X3_FAR60921_0ea3ede7/fastq_pass
### Set the path to look in for failed fastq files
failpath=/home/uqsidris/network/SEQLAB-Q4547/cDNA-PCR_mRNA_vaccine_QC_010/mRNA_vaccines/20220422_0640_X3_FAR60921_0ea3ede7/fastq_fail
### Set the path to save outputs
outputpath=/home/uqsidris/workdir/ReanalysisVac14_cDNA_20220422_0640_X3_FAR60921_0ea3ede7
### Set the path to look in for the reference fasta files
referencepath=/home/uqsidris/workdir/referencefiles
### Set path to seq run summary file
runsummaryfilepath=/home/uqsidris/network/SEQLAB-Q4547/cDNA-PCR_mRNA_vaccine_QC_010/mRNA_vaccines/20220422_0640_X3_FAR60921_0ea3ede7/sequencing_summary_FAR60921_24a6bd1f
### Set the path to look in for the bed files
bedpath=/home/uqsidris
### Set array for the s used match format to the name of the folders that the data is in
barcode_arr=(barcode08)
### Set array for the sample names that match in the order of the s in the array above
sample_arr=(2022_14_482162_F4_covid_vac_base)
### Set array for the reference files in fasta format that match in the order of the samples in the array above
reference_arr=(2022_14_482162_F4)
### Set array for the bed files (mRNA start coordinate eg cap site, kozak or orf start) that match in the order of the samples in the array above
startbed_arr=(14_kozak_start)
### Set array for the bed files (end coordinate 3utr) that match in the order of the samples in the array above
endbed_arr=(14_3utr_end)

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
    zcat $passpath/$barcodefile/*fastq.gz > $outputpath/$samplefile/$samplefile'_allpassedreads.fastq.gz'
done

# Join failed fastq files from each barcode then rename and save in appropriate sample folder
for i in ${!barcode_arr[@]}; do
barcodefile=${barcode_arr[i]}
samplefile=${sample_arr[i]}
zcat $failpath/$barcodefile/*fastq.gz > $outputpath/$samplefile/$samplefile'_allfailedreads.fastq.gz'
done

# change file name to .fastq
for i in ${!barcode_arr[@]}; do
    barcodefile=${barcode_arr[i]}
    samplefile=${sample_arr[i]}
    cat $outputpath/$samplefile/$samplefile'_allpassedreads.fastq.gz' > $outputpath/$samplefile/$samplefile'_allpassedreads.fastq'
    cat $outputpath/$samplefile/$samplefile'_allfailedreads.fastq.gz' > $outputpath/$samplefile/$samplefile'_allfailedreads.fastq'
done

#pychopper
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    cdna_classifier.py -r $outputpath/$samplefile/$samplefile'_passedreads_pychop_oldcdn.pdf' -u $outputpath/$samplefile/$samplefile'_passedreads_pychop_oldcdn_unclassified.fastq' -w $outputpath/$samplefile/$samplefile'_passedreads_pychop_oldcdn_rescued.fastq' -S $outputpath/$samplefile/$samplefile'_passedreads_pychop_oldcdn_stats.tsv' $outputpath/$samplefile/$samplefile'_allpassedreads.fastq' $outputpath/$samplefile/$samplefile'_allpassedreads_pychop_oldcdn_full_length.fastq'
    cdna_classifier.py -r $outputpath/$samplefile/$samplefile'_failedreads_pychop_oldcdn.pdf' -u $outputpath/$samplefile/$samplefile'_failedreads_pychop_oldcdn_unclassified.fastq' -w $outputpath/$samplefile/$samplefile'_failedreads_pychop_oldcdn_rescued.fastq' -S $outputpath/$samplefile/$samplefile'_failedreads_pychop_oldcdn_stats.tsv' $outputpath/$samplefile/$samplefile'_allfailedreads.fastq' $outputpath/$samplefile/$samplefile'_allfailedreads_pychop_oldcdn_full_length.fastq'
done

# Join rescued and full length fastq files made by pychopper from each barcode then rename and save in appropriate sample folder
for i in ${!barcode_arr[@]}; do
    barcodefile=${barcode_arr[i]}
    samplefile=${sample_arr[i]}
    cat $outputpath/$samplefile/$samplefile'_passedreads_pychop_oldcdn_rescued.fastq' $outputpath/$samplefile/$samplefile'_allpassedreads_pychop_oldcdn_full_length.fastq' > $outputpath/$samplefile/$samplefile'_allpassedreads_pychop_oldcdn_joined.fastq'
    cat $outputpath/$samplefile/$samplefile'_failedreads_pychop_oldcdn_rescued.fastq' $outputpath/$samplefile/$samplefile'_allfailedreads_pychop_oldcdn_full_length.fastq' > $outputpath/$samplefile/$samplefile'_allfailedreads_pychop_oldcdn_joined.fastq'
done

# Running minimap2 align fastq to reference fasta and make sam file in correct directory for all passed reads
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    minimap2 -ax map-ont $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_pychop_oldcdn_joined.fastq' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn.sam'
    minimap2 -ax map-ont $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allfailedreads_pychop_oldcdn_joined.fastq' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn.sam'
done

-----------------------------------------------------------
                    MANA STARTS HERE
-----------------------------------------------------------

# Run samtools for various statistics and also create a sorted and indexed bam file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # create samtools stats file for the sam file for all passed reads
    samtools stats $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn.sam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_stats.txt'
    samtools stats $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn.sam' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_stats.txt'
    # create a bam file from the sam file for all passed reads
    samtools view -S -b $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn.sam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn.bam'
    samtools view -S -b $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn.sam' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn.bam'
    # sort the bam file for all passed reads
    samtools sort $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam'
    samtools sort $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn.bam' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted.bam'
    # index the bam file for all passed reads
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam'
    samtools index $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted.bam'
    # create samtools flagstat file for all passed reads
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_flagstat.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted.bam' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted_flagstat.txt'
done

# Run pysamstats to get various errors per base then use awk to remove obsolete columns and create new file then do some calculations and output a new file
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    # make pysamstats txt file
    pysamstats --max-depth=300000000 --fasta $referencepath/$referencefile'.fasta' --type variation $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_pysam.txt'
    # all together cut and percent each error in new column with headers
    cut -f 5,7,9,11,13,15,17,19,21,23 --complement $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR==1{$13 = $13 OFS "Match per base" OFS "Mismatch per base" OFS "Deletion per base" OFS "Insertion per base"} NR>1{$13 = $13 OFS ($5/$4) OFS ($6/$4) OFS ($7/$4) OFS ($8/$4)}1' CONVFMT='%.3f' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt'
    # making a table with min max average and headers.
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$14}; if($14>max) {max=$14}; if($14<min) {min=$14}; total+=$14; count+=1} {CONVFMT="%.3f"} END {print "Error Performance per Nucleotide" "\n" "Frequency %" OFS "avg" OFS "max" OFS "min" "\n" "Match" OFS total/count OFS max OFS min}' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt'
    # Next Line next error insertion
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$17}; if($17>max) {max=$17}; if($17<min) {min=$17}; total+=$17; count+=1} {CONVFMT="%.3f"} END {print "Insertion" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt'
    # Next Line next error mismatch
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$15}; if($15>max) {max=$15}; if($15<min) {min=$15}; total+=$15; count+=1} {CONVFMT="%.3f"} END {print "Mismatch" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt'
    # Next Line next error deletion
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$16}; if($16>max) {max=$16}; if($16<min) {min=$16}; total+=$16; count+=1} {CONVFMT="%.3f"} END {print "Deletion" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt'
    # Next Line next total error
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR>4{if(min==""){min=$4}; if(max==""){max=$3}; if ($3>max) {max=$3}; if($4<min) {min=$4}; total+=$2; count+=1} {CONVFMT="%.3f"} END {print "Total Error" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt'
    # Notes
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt' | awk 'BEGIN{FS=OFS="\t"} END {print "Notes" "\n" "Error = mismatch + deletion" "\n" "Mismatch does not include deletion"}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_pysamreporttable.txt'
done

# Run awk etc to get various numbers from flagstat, stats and pysamtable do some calculations and output a new file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # making a table with title Run performance
    echo "Run Performance" > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    # convert flagstat file to tab delimited file
    cat $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted_flagstat.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted_flagstat_tabdelimited.txt'
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_flagstat.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_flagstat_tabdelimited.txt'
    # open tab delimited flagstat file and print total reads number (1st number in the file) for passed to the report
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$1} END {print "Passed Reads >Q9" OFS passedread}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{failedread=$1} END {print "Failed Reads <Q9" OFS failedread}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    # sum the above to get total reads from the run
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} NR==3{failedread=$2} END {print "Total Reads" OFS passedread+failedread}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    # find flagstat mapped number 7th number in file
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==8{mappedread=$1} END {print "Mapped Reads" OFS mappedread}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    # calc unmapped reads 
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} NR==5{mappedread=$2} END {print "Unmapped Reads" OFS passedread-mappedread}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    # convert samtools stats file to tab delimited file
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_stats.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_stats_tabdelimited.txt'
    # find read length average from samtools stats
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_stats_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==33{readlengthavg=$4} END {print "Read Length Average" OFS readlengthavg}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    # find coverage -average -max -min from the pysamstats report
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Average Coverage" OFS total/count}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Min Coverage" OFS min}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Max Coverage" OFS max}' >> $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_allreads_runperformance_reporttable.txt'
done

# Create file listing the lengths of each primary mapped read in the sorted bam file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools view -F 2048 $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam' | awk '{print length($10)}'| sort -n | uniq -c |awk ' { t = $1; $1 = $2; $2 = t; print; } '|tr ' ' '\t'|sed '1d' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_primarymappedreadlengths.txt'
    samtools view -F 2048 $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted.bam' | awk '{print length($10)}'| sort -n | uniq -c |awk ' { t = $1; $1 = $2; $2 = t; print; } '|tr ' ' '\t'|sed '1d' > $outputpath/$samplefile/$samplefile'_failed_pychop_oldcdn_sorted_primarymappedreadlengths.txt'
done

# Check full length and degraded from 3 or 5 prime ends
for i in ${!startbed_arr[@]}; do
    startbedfile=${startbed_arr[i]}
    samplefile=${sample_arr[i]}
    ##### intersect all reads bam with reads overlapping at start coordinate from start.bed file make a bam and index and flagstat
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam' -b $bedpath/$startbedfile'.bed' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_startcoord.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_startcoord.bam'
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_startcoord.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_startcoord_flagstat.txt'
done



##### Get full length bam --intersect newly made start coord bam with reads ending at last base 3utr coordinate from endutr.bed file ---- index and flagstat
for i in ${!endbed_arr[@]}; do
    endbedfile=${endbed_arr[i]}
    samplefile=${sample_arr[i]}
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_startcoord.bam' -b $bedpath/$endbedfile'.bed' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_full_length_reads.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_full_length_reads.bam'
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_full_length_reads.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_full_length_reads_flagstat.txt'

    ##### Get possibly degraded has a start but No_3utr_overlap reads bam --intersect start coord bam with endutr.bed file but pull out the reads that DONT INTERSECT option -v  ---- index and flagstat
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_startcoord.bam' -b $bedpath/$endbedfile'.bed' -v > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_start_No_3utr_overlap_reads.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_start_No_3utr_overlap_reads.bam'
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_start_No_3utr_overlap_reads.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_start_No_3utr_overlap_reads_flagstat.txt'

    ##### intersect all reads bam with reads overlapping at 3utr end coordinate from end.bed file make a bam and index and flagstat
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam' -b $bedpath/$endbedfile'.bed' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_3utrendcoord.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_3utrendcoord.bam'
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_3utrendcoord.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_3utrendcoord_flagstat.txt'
done

##### Get possibly degraded No_start_overlap but has a end 3utr reads bam --intersect start coord bam with endutr.bed file but pull out the reads that DONT INTERSECT option -v  ---- index and flagstat
for i in ${!startbed_arr[@]}; do
    startbedfile=${startbed_arr[i]}
    samplefile=${sample_arr[i]}
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_3utrendcoord.bam' -b $bedpath/$startbedfile'.bed' -v > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_3utr_No_start_overlap_reads.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_3utr_No_start_overlap_reads.bam'
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_3utr_No_start_overlap_reads.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_3utr_No_start_overlap_reads_flagstat.txt'
done

# Check off target
for i in ${!startbed_arr[@]}; do
    startbedfile=${startbed_arr[i]}
    samplefile=${sample_arr[i]}
    ##### intersect all reads bam with reads NOT overlapping at start coordinate from start.bed file make a bam and index and flagstat
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted.bam' -b $bedpath/$startbedfile'.bed' -v > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord.bam'
    samtools flagstat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord_flagstat.txt'
done

##### Get reads NOT overlapping at start coordinate bam --intersect with last base 3utr coordinate from endutr.bed file Get reads NO overlap at 3utr ---- index and flagstat ---this should be completely off target reads.
for i in ${!endbed_arr[@]}; do
    endbedfile=${endbed_arr[i]}
    samplefile=${sample_arr[i]}
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord.bam' -b $bedpath/$endbedfile'.bed' -v > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord_NO_3utr_end_coord.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord_NO_3utr_end_coord.bam'
    samtools flagstat  $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord_NO_3utr_end_coord.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord_NO_3utr_end_coord_flagstat.txt'

# pulling out the reads
###total_mapped_reads all reads
cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_sorted_flagstat.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i Total_reads\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_TotalReads.txt'
###full length_reads_start and stop overlap correct reads
cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_full_length_reads_flagstat.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i full_length_reads\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_FullLength.txt'
###no start but has 3utr end overlap degraded
cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_3utr_No_start_overlap_reads_flagstat.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i degraded_from5prime\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_degraded_from5prime.txt'
###has start but no 3utr end overlap degraded
cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_Has_start_No_3utr_overlap_reads_flagstat.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i degraded_from3prime\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_degraded_from3prime.txt'
###no start or 3utr end overlap so completely offtarget or other degraded types
cat $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_NO_startcoord_NO_3utr_end_coord_flagstat.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i off_target_other_degraded_reads\t' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_offtarget_otherdegraded.txt'
# Compiling all reads into one file
paste $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_TotalReads.txt' $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_FullLength.txt' $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_degraded_from5prime.txt' $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_degraded_from3prime.txt' $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_offtarget_otherdegraded.txt' | column -t > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_reported_reads.txt'
done

# Make a bam of unmapped reads
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools view -S -b -f 4 $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn.sam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_unmapped_passedreads.bam'
    samtools index $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_unmapped_passedreads.bam'
done

# Converting unmapped bam to unmapped fasta to be able to use blastn website?
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools fasta $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_unmapped_passedreads.bam' > $outputpath/$samplefile/$samplefile'_passed_pychop_oldcdn_unmapped_passedreads.fasta'
done

#Add variant calling and consensus generation
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' | bcftools call -c -M --ploidy 1 -Oz -o $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz'
    bcftools index $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz'
    bcftools norm -f $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz' -Ob -o $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_BCF.bcf'
    bcftools consensus -a - -f $referencepath/$referencefile'.fasta' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_VCF.vcf.gz' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_BCF_consensus.fa'
done

#!/bin/bash
echo "Done"
