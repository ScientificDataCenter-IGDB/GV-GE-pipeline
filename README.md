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
7.	VCFtools (v0.1.15)
8.	vcflib:Vcffilter (v1.0.9)
9.	bcftools (v1.9)
10.	annovar(version 2020-06-08)
11.	plink (v1.90b6.21) 
12.	bedtools(v2.30.0)
13.	python(v3.7.6)
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
## For trio:
5. Calling
```
python Platypus.py callVariants --bamFiles=off_sorted_dedup.bam,mo_sorted_dedup.bam,fa_sorted_dedup.bam --refFile=ref.fa  --output=trio.PlatypusTrio.vcf --maxVariants=20 --maxReads=30000000 --nCPU=20 --logFileName=Trio.log.txt
```
6. Select SNPs and INDELs
```
gatk --java-options -Xmx20G SelectVariants -V trio.PlatypusTrio.vcf -select-type SNP -O trio.SNP.PlatypusTrio.vcf
gatk --java-options -Xmx20G SelectVariants -V trio.PlatypusTrio.vcf -select-type INDEL -O trio.INDEL.PlatypusTrio.vcf
```
7. Site filtering
```
