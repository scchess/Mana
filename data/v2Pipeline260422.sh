# Test bash shell for complete binf pipeline

# Set the path where this bash shell will be
shellpath=/home/uqsidris/'v2Pipeline260422.sh'

# Set the path to look in for passed fastq files
passpath=/home/uqsidris/DATA_RUNS/cDNA-PCR_mRNA_vaccine_QC_010/mRNA_vaccines/20220422_0640_X3_FAR60921_0ea3ede7/fastq_pass

# Set the path to look in for failed fastq files
failpath=/home/uqsidris/DATA_RUNS/cDNA-PCR_mRNA_vaccine_QC_010/mRNA_vaccines/20220422_0640_X3_FAR60921_0ea3ede7/fastq_fail

# Set the path to save outputs
outputpath=/home/uqsidris/workdir/20220422_0640_X3_FAR60921_0ea3ede7

# Set the path to look in for the reference fasta files
referencepath=/home/uqsidris/workdir/cDNA007

# Set path to seq run summary file
runsummaryfilepath=/home/uqsidris/DATA_RUNS/cDNA-PCR_mRNA_vaccine_QC_010/mRNA_vaccines/20220422_0640_X3_FAR60921_0ea3ede7/sequencing_summary_FAR60921_24a6bd1f

# Set the path to look in for the bed files
bedpath=/home/uqsidris

# Set array for the Barcodes used match format to the name of the folders that the data is in
barcode_arr=(barcode11 barcode12)

# Set array for the sample names that match in the order of the barcodes in the array above
sample_arr=(cDNA_UnMod_37C_NEBT7_BaseGfpmRNA_polyA cDNA_Mod_37C_NEBT7_BaseGfpmRNA_polyA)

# Set array for the reference files in fasta format that match in the order of the samples in the array above
reference_arr=(plasmid_gfp_ref plasmid_gfp_ref)

# Set array for the bed files that match in the order of the samples in the array above
bed_arr=(Target_PolyA_Overlap Target_PolyA_Overlap)

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

# Join passed fastq files from each barcode then rename and save in appropriate sample folder
for i in ${!barcode_arr[@]}; do
barcodefile=${barcode_arr[i]}
samplefile=${sample_arr[i]}
zcat $passpath/$barcodefile/*fastq.gz > $outputpath/$samplefile/$samplefile'_allpassedreads.fastq.gz'
done

# FastQC
for i in ${!sample_arr[@]}; do
samplefile=${sample_arr[i]}
cat $outputpath/$samplefile/$samplefile'_allpassedreads.fastq.gz' > $outputpath/$samplefile/$samplefile'_allpassedreads.fastq'
fastqc -t 14 -o $outputpath/$samplefile/ $outputpath/$samplefile/$samplefile'_allpassedreads.fastq'
done

# Join failed fastq files from each barcode then rename and save in appropriate sample folder
for i in ${!barcode_arr[@]}; do
barcodefile=${barcode_arr[i]}
samplefile=${sample_arr[i]}
zcat $failpath/$barcodefile/*fastq.gz > $outputpath/$samplefile/$samplefile'_allfailedreads.fastq.gz'
done

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
    # create samtools coverage file for all passed reads
    samtools coverage $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_coverage.txt'
    samtools coverage $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_coverage.txt'
    # create samtools flagstat file for all passed reads
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat.txt'
    # create samtools depth file for all passed reads
    samtools depth -a $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_depth.txt'
    samtools depth -a $outputpath/$samplefile/$samplefile'_allfailedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_depth.txt'
done

# Run pysamstats to get various errors per base then use awk to remove obsolete columns and create new file then do some calculations and output a new file
for i in ${!reference_arr[@]}; do
    referencefile=${reference_arr[i]}
    samplefile=${sample_arr[i]}
    # make pysamstats txt file
    pysamstats --max-depth=3000000 --fasta $referencepath/$referencefile'.fasta' --type variation $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_pysam.txt'
    # all together cut and percent each error in new column with headers
    cut -f 5,7,9,11,13,15,17,19,21,23 --complement $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR==1{$13 = $13 OFS "Match per base" OFS "Mismatch per base" OFS "Deletion per base" OFS "Insertion per base"} NR>1{$13 = $13 OFS ($5/$4) OFS ($6/$4) OFS ($7/$4) OFS ($8/$4)}1' CONVFMT='%.3f' > $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt'
    # making a table with min max average and headers.
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$14}; if($14>max) {max=$14}; if($14<min) {min=$14}; total+=$14; count+=1} {CONVFMT="%.3f"} END {print "Error Performance per Nucleotide" "\n" "Frequency %" OFS "avg" OFS "max" OFS "min" "\n" "Match" OFS total/count OFS max OFS min}' > $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt'
    # Next Line next error mismatch
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$15}; if($15>max) {max=$15}; if($15<min) {min=$15}; total+=$15; count+=1} {CONVFMT="%.3f"} END {print "Mismatch" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt'
    # Next Line next error deletion
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$16}; if($16>max) {max=$16}; if($16<min) {min=$16}; total+=$16; count+=1} {CONVFMT="%.3f"} END {print "Deletion" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt'
    # Next Line next error insertion
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$17}; if($17>max) {max=$17}; if($17<min) {min=$17}; total+=$17; count+=1} {CONVFMT="%.3f"} END {print "Insertion" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt'
    # Next Line next total error
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR>3{if(min==""){min=$4}; if(max==""){max=$3}; if ($3>max) {max=$3}; if($4<min) {min=$4}; total+=$2; count+=1} {CONVFMT="%.3f"} END {print "Total Error" OFS total/count OFS max OFS min}' >> $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt'
    # Notes
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt' | awk 'BEGIN{FS=OFS="\t"} END {print "Notes" "\n" "Error = mismatch + insertion + deletion" "\n" "Mismatch does not include deletion"}' >> $outputpath/$samplefile/$samplefile'_allpassedreads_pysamreporttable.txt'
done

# Run awk etc to get various numbers from flagstat, stats and pysamtable do some calculations and output a new file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # making a table with title Run performance
    echo "Run Performance" > $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    # convert flagstat file to tab delimited file
    cat $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat_tabdelimited.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat_tabdelimited.txt'
    # open tab delimited flagstat file and print total reads number (1st number in the file) for passed to the report
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==1{passedread=$1} END {print "Passed Reads >Q9" OFS passedread}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_allfailedreads_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==1{failedread=$1} END {print "Failed Reads <Q9" OFS failedread}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    # sum the above to get total reads from the run
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} NR==3{failedread=$2} END {print "Total Reads" OFS passedread+failedread}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    # find flagstat mapped number 7th number in file
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_flagstat_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==7{mappedread=$1} END {print "Mapped Reads" OFS mappedread}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    # calc unmapped reads 
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} NR==5{mappedread=$2} END {print "Unmapped Reads" OFS passedread-mappedread}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    
    # convert samtools stats file to tab delimited file
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_stats.txt' | tr ' ' '\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_stats_tabdelimited.txt'
    # find read length average from samtools stats
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_stats_tabdelimited.txt' | awk 'BEGIN{FS=OFS="\t"} NR==33{readlengthavg=$4} END {print "Read Length Average" OFS readlengthavg}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    # find coverage -average -max -min from the pysamstats report
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Average Coverage" OFS total/count}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Min Coverage" OFS min}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_edited_pysam.txt' | awk 'BEGIN{FS=OFS="\t"} NR>1{if(min==""){min=max=$4}; if($4>max) {max=$4}; if($4<min) {min=$4}; total+=$4; count+=1} {CONVFMT="%.3f"} END {print "Max Coverage" OFS max}' >> $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt'
done

# Make a table with fractions for each sample to subsample to varying levels so can enter these into the next shell
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # making a table with title Run performance
    echo "Subsampling fractions to hit certain read targets" > $outputpath/$samplefile/$samplefile'_subsample_fraction_calc_table.txt'
    # look at number of passed reads in the report and calculate what is fraction for 200,000 read subsample
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} END {print "fraction to get 200k reads" OFS 200000/passedread}' >> $outputpath/$samplefile/$samplefile'_subsample_fraction_calc_table.txt'
    # look at number of passed reads in the report and calculate what is fraction for 100,000 read subsample
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} END {print "fraction to get 100k reads" OFS 100000/passedread}' >> $outputpath/$samplefile/$samplefile'_subsample_fraction_calc_table.txt'
    # look at number of passed reads in the report and calculate what is fraction for 50,000 read subsample
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} END {print "fraction to get 50k reads" OFS 50000/passedread}' >> $outputpath/$samplefile/$samplefile'_subsample_fraction_calc_table.txt'
    # look at number of passed reads in the report and calculate what is fraction for 10,000 read subsample
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} END {print "fraction to get 10k reads" OFS 10000/passedread}' >> $outputpath/$samplefile/$samplefile'_subsample_fraction_calc_table.txt'
    # look at number of passed reads in the report and calculate what is fraction for 5,000 read subsample
    cat $outputpath/$samplefile/$samplefile'_allreads_runperformance_reporttable.txt' | awk 'BEGIN{FS=OFS="\t"} NR==2{passedread=$2} END {print "fraction to get 5k reads" OFS 5000/passedread}' >> $outputpath/$samplefile/$samplefile'_subsample_fraction_calc_table.txt'
done

# Create file listing the lengths of each mapped read in the sorted bam file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools view -F 4 $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' | awk '{print length($10)}'| sort -n | uniq -c |awk ' { t = $1; $1 = $2; $2 = t; print; } '|tr ' ' '\t'|sed '1d' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_mappedreadlengths.txt'
done

# Check on and off target
for i in ${!bed_arr[@]}; do
    bedfile=${bed_arr[i]}
    samplefile=${sample_arr[i]}
    # on target reads
    bedtools intersect -a $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' -b $bedpath/$bedfile'.bed' -f 0.98 -r > $outputpath/$samplefile/$samplefile'_allpassedreads_ontarget.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_ontarget.bam'
    # off target reads
    bedtools intersect -abam $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' -b $outputpath/$samplefile/$samplefile'_allpassedreads_ontarget.bam' -v -f 1.0 -r > $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.bam'
    # filtering off target reads to different categories
    samtools view -b $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.bam' "2021-8:1-431" > $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_upstream.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_upstream.bam'
    samtools view -b $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.bam' "2021-8:1636-3790" > $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_downstream.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_downstream.bam'
    bedtools intersect -a off_target.bam -b $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_upstream.bam' -b $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_downstream.bam' -v -f 1.0 -r > $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_truncated.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_truncated.bam'
    # flagstat above files
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_ontarget.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_ontarget.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_downstream.bam' > $$outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_unlinearised.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_upstream.bam' >> $$outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_unlinearised.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_truncated.bam' > $$outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_truncated.txt'
    samtools flagstat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' > $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_total_mapped_reads.txt'
    # pulling out the reads
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_sorted_total_mapped_reads.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i Total_reads\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_full.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_ontarget.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i on_target_reads\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_on.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i off_target_reads\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_off.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_unlinearised.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} NR==24{primarymappedreads2=$1} END {print primarymappedreads+primarymappedreads2}'| sed '1i off_target_unlinearised_reads\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_off_uplinearised.txt'
    cat $outputpath/$samplefile/$samplefile'_allpassedreads_offtarget_truncated.txt' |tr ' ' '\t'| awk 'BEGIN{FS=OFS="\t"} NR==8{primarymappedreads=$1} END {print primarymappedreads}'| sed '1i off_target_truncated_reads\t' > $outputpath/$samplefile/$samplefile'_allpassedreads_off_truncated.txt'
    # Compiling all reads into one file
    paste $outputpath/$samplefile/$samplefile'_allpassedreads_full.txt' $outputpath/$samplefile/$samplefile'_allpassedreads_on.txt' $outputpath/$samplefile/$samplefile'_allpassedreads_off.txt' $outputpath/$samplefile/$samplefile'_allpassedreads_off_uplinearised.txt' $outputpath/$samplefile/$samplefile'_allpassedreads_off_truncated.txt' | column -t > $outputpath/$samplefile/$samplefile'_allpassedreads_reported_reads.txt'
done

# Make a bam of unmapped reads
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools view -S -b -f 4 $outputpath/$samplefile/$samplefile'_allpassedreads.sam' > $outputpath/$samplefile/$samplefile'_unmapped_passedreads.bam'
    samtools index $outputpath/$samplefile/$samplefile'_unmapped_passedreads.bam'
done

# Converting unmapped bam to unmapped fasta to be able to use blastn website?
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    samtools fasta $outputpath/$samplefile/$samplefile'_unmapped_passedreads.bam' > $outputpath/$samplefile/$samplefile'_unmapped_passedreads.fasta'
    #cat ~/workdir/20220329_2343_X4_FAR61095_281da2bc/cDNA_Mod_32C_NEBT7_BaseGfpmRNA/cDNA_Mod_32C_NEBT7_BaseGfpmRNA_unmapped_passedreads.bam | awk '{printf(">%s\n%s\n",$1,$10);}' | blastn -db contamination -out ~/workdir/20220329_2343_X4_FAR61095_281da2bc/cDNA_Mod_32C_NEBT7_BaseGfpmRNA/cDNA_Mod_32C_NEBT7_BaseGfpmRNA_unmapped_passedreads.blast
done

# Separating forward and reverse reads to new bams then counting them and making result output file
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    # Using sam flag to count forward reads from the original bam and make a forward reads only bam
    samtools view -F 16 -o $outputpath/$samplefile/$samplefile'_allpassedreads_onlyForward.bam' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam'
    # Using sam flag to count reverse reads from the original bam and make a reverse reads only bam
    samtools view -f 16 -o $outputpath/$samplefile/$samplefile'_allpassedreads_onlyReverse.bam' $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam'
    samtools index $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam'
    # then count words in bam to count forward reads and make a new file
    samtools view $outputpath/$samplefile/$samplefile'_allpassedreads_onlyForward.bam' |wc -l > $outputpath/$samplefile/$samplefile'_allpassedreads_ForwardReadsCount.txt'
    # then count words in bam to count forward reads and make a new file
    samtools view $outputpath/$samplefile/$samplefile'_allpassedreads_onlyReverse.bam' |wc -l > $outputpath/$samplefile/$samplefile'_allpassedreads_ReverseReadsCount.txt'
done

# nanoplot for the whole run and each sample
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    NanoPlot --summary $runsummaryfilepath'.txt' --loglength -o $outputpath/$samplefile/$samplefile'run_nanoplot'
    NanoPlot --color green --bam $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' -o $outputpath/$samplefile/$samplefile'sample_nanoplot'
done

# pycoQC for the whole run and each sample
for i in ${!sample_arr[@]}; do
    samplefile=${sample_arr[i]}
    pycoQC -f $runsummaryfilepath'.txt' -a $outputpath/$samplefile/$samplefile'_allpassedreads_sorted.bam' -o $outputpath/$samplefile/$samplefile'run_pycoQC.html'
done

#run most of the above scripts again for subsampling put in a different shell as not sure how to integrate here...need to figure out best read depth to subsample to....
#add variant calling and consensus fasta generator.
#look into short fragments on poly A section
#look into the cascade from larger to shorter reads across the length of the reference, something to do with reverse transcription or strand switching? Compare to dRNA
#pyycotorch cdna use for the input to mapping? as it makes a fastq set with just reads containing both adaptors and trims the adaptors orientates everything to forward reads..
#or figure out our own way to manually grep and filter the fastqs for the adaptors
#minimap with alternate setting for cdna
#using only primary mapped reads in igv?
#check unmapped for adaptors - map to adaptor sequences
#binning lenths graph - still to do
#blast is installed but no database..not sure how to get that working
#figure out how to make a log file of versions of tools used in the pipeline
#Swatis on and off target section might need to be further refined for dRNA
