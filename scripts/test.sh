
sudo docker build -t mana .
sudo docker run -i -t -v ${PWD}:/data mana /bin/bash
cd /data
nanosim-h GCF_000005845.2_ASM584v2_genomic.fna
mv simulated.fa simulated.fastq
minimap2 -ax map-ont GCF_000005845.2_ASM584v2_genomic.fna simulated.fastq | samtools view -bS | samtools sort > simulated.bam
samtools view simulated.bam | wc
samtools index simulated.bam
samtools flagstat simulated.bam
pysamstats --max-depth=300000000 --fasta GCF_000005845.2_ASM584v2_genomic.fna --type variation simulated.bam > simulated.bam.txt

bcftools mpileup -d 300000000 --no-BAQ --min-BQ 0 -Ou -f GCF_000005845.2_ASM584v2_genomic.fna simulated.bam | bcftools call -c -M --ploidy 1 -Oz -o simulated.vcf.gz
bcftools index simulated.vcf.gz
bcftools norm -f GCF_000005845.2_ASM584v2_genomic.fna $simulated.vcf.gz -Ob -o simulated.bcf
bcftools consensus -a - -f GCF_000005845.2_ASM584v2_genomic.fna simulated.vcf.gz > simulated.consensus.fa
