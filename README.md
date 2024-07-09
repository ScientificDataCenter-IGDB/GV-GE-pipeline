# GV-GE-pipeline
# 1. Introduction
--------------------
This pipeline aims to detect de novo mutations(DNMs) in the pedigrees of large animals with gene-edited, filter and annotate DNMs accordingly.
# 2.Software version
--------------------
1.  fastp (v0.23.1)
2.	bwa-mem2 (v2.2.1)
3.	Qulimap(v2.2.1)
4.	samtools (v1.19.2)
5.	GATK (v4.2.6.1)
6.	Platypus (v0.8.1)
7.	vcflib:Vcffilter (v1.0.9)
8.	annovar(version 2020-06-08)
9.	python(v2.7.5)
10.	perl(v5.34.0)
# 3.Usage
--------------------
## For single sample:
1. Processing fastq files
```
fastp -i R1.fq.gz -I R2.fq.gz -o R1_qc.fq.gz -O R2_qc.fq.gz --html SAMPLE.html --json SAMPLE.json
```
2. Mapping
```
bwa-mem2 mem -M -t 4 -R "@RG\tID:SAMPLEtSM:SAMPLE\tLB:SAMPLE\tPL:ILLUMINA\tPU:run"  ref.fa R1_qc.fq.gz R2_qc.fq.gz > SAMPLE.sam 2> SAMPLE.sam.err
```
3. sam->bam & sort
```
gatk --java-options -Xmx20G SortSam -I SAMPLE.sam -O SAMPLE_sorted.bam -SO coordinate --CREATE_INDEX true --VALIDATION_STRINGENCY LENIENT
```
4. Mark Duplicates
```
gatk --java-options -Xmx20G MarkDuplicates -I SAMPLE_sorted.bam -O SAMPLE_sorted_dedup.bam -M SAMPLE_metrics.txt --ASSUME_SORTED true --VALIDATION_STRINGENCY LENIENT --CREATE_INDEX true
```
5. Depth caculation
```
qualimap bamqc -bam SAMPLE_sorted_dedup.bam -outdir ./ -outformat PDF:HTML --java-mem-size=50G --skip-duplicated
```
## For trio:
6. Calling
```
python Platypus.py callVariants --bamFiles=off_sorted_dedup.bam,mo_sorted_dedup.bam,fa_sorted_dedup.bam --refFile=ref.fa  --output=trio.PlatypusTrio.vcf --maxVariants=20 --maxReads=30000000 --nCPU=20 --logFileName=Trio.log.txt
```
7. Separate SNPs and INDELs
```
gatk --java-options -Xmx20G SelectVariants -V trio.PlatypusTrio.vcf -select-type SNP -O trio.SNP.PlatypusTrio.vcf
gatk --java-options -Xmx20G SelectVariants -V trio.PlatypusTrio.vcf -select-type INDEL -O trio.INDEL.PlatypusTrio.vcf
```
8. Site filtering
```
vcffilter -f "SbPval >0.001 & HapScore <13 & MQ >30 & QD >2 & NR >0 & NF >0" trio.SNP.PlatypusTrio.vcf > trio.SNP.PlatypusTrio.filterSF.vcf
vcffilter -f "SbPval >0.001 & HapScore <13 & MQ >30 & QD >2 & NR >0 & NF >0" trio.INDEL.PlatypusTrio.INDEL.vcf > trio.INDEL.PlatypusTrio.filterSF.vcf
```
9. Genotype and vaf filtering
```
avg_depth=$(grep "mean coverageData = " ./qualimap/genome_results.txt  |perl -pe 's/\s+mean coverageData = //' |perl -pe 's/X//')
low=$(($avg_depth/3))
dp_low=$(($low+1))
dp_high=$(($avg_depth*2))
perl filter_genotype.pl --input trio.SNP.PlatypusTrio.filterSF.vcf --dp $dp_low --dp_up $dp_high -gq 40 --out trio.SNP.PlatypusTrio.filterSF_GF.vcf
perl filter_genotype.pl --input trio.INDEL.PlatypusTrio.filterSF.vcf --dp 5 --dp_up $dp_high -gq 5 --out trio.INDEL.PlatypusTrio.filterSF_GF.vcf
```
10. defined DNMs
```
perl defined_DNMs.pl --input trio.SNP.PlatypusTrio.filterSF_GF.vcf --samplename $off_sample_name --vaf 0.2 --out trio.SNP.DNMs.vcf
perl defined_DNMs.pl --input trio.INDEL.PlatypusTrio.filterSF_GF.vcf --samplename $off_sample_name --vaf 0.2 --out trio.INDEL.DNMs.vcf
```
11. DNMs filtering
```
grep -P "#|\tPASS\t" trio.SNP.DNMs.vcf > trio.SNP.DNMs.PASS.vcf
perl DNMs_SNP_filtering.pl --input trio.SNP.DNMs.PASS.vcf --out trio.SNP.DNMs.PASS.filtered.vcf --off $off_sample_name --fa $fa_sample_name --mo $mo_sample_name --spe ref.fa
grep -P "#|\tPASS\t" trio.INDEL.DNMs.vcf > trio.INDEL.DNMs.PASS.vcf
perl DNMs_INDEL_filtering.pl --input trio.INDEL.DNMs.PASS.vcf --out trio.INDEL.DNMs.PASS.filtered.vcf --off $off_sample_name --fa $fa_sample_name --mo $mo_sample_name --spe ref.fa
```
12. DNMs annotation
```
perl /PATH/annovar/table_annovar.pl trio.SNP.DNMs.PASS.filtered.vcf /PATH/SPECIES/SPECIESdb -buildver SPECIES  -out trio.SNP.DNMs.PASS.filtered -remove -protocol refGene -operation g -nastring . --nopolish -vcfinput
perl /PATH/annovar/table_annovar.pl trio.INDEL.DNMs.PASS.filtered.vcf /PATH/SPECIES/SPECIESdb -buildver SPECIES  -out trio.INDEL.DNMs.PASS.filtered -remove -protocol refGene -operation g -nastring . --nopolish -vcfinput
```
