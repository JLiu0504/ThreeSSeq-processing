#!/bin/bash

ref=NC_000913.3.fa
sample=776input-2_S9
echo $sample

bwa mem ${ref} ${sample}_R1.trim.fastq ${sample}_L001_R2_001.fastq | samtools view -bhS - > ${sample}_aln.map.bam
samtools sort ${sample}_aln.map.bam -o ${sample}_aln.sort.bam
samtools index ${sample}_aln.sort.bam

java -jar /home/jingjing/Downloads/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${sample}_aln.sort.bam O=${sample}_aln.rmdup.bam M=${sample}_dup_metrics.txt 
samtools index ${sample}_aln.rmdup.bam

