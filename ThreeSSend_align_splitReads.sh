#!/bin/bash

ref=NC_000913.3.fa
sample=ExID-FLAG_rpoB2_2
trimR1=${sample}_R1.trim.fastq

gunzip ${sample}_R1_001.fastq.gz
gunzip ${sample}_R2_001.fastq.gz

java -jar /home/jingjing/Downloads/Trimmomatic-0.36/trimmomatic-0.36.jar SE ${sample}_R1_001.fastq $trimR1 HEADCROP:3

bwa mem ${ref} $trimR1 ${sample}_R2_001.fastq | samtools view -bhS - > ${sample}_aln.map.bam
samtools sort ${sample}_aln.map.bam -o ${sample}_aln.sort.bam
samtools index ${sample}_aln.sort.bam

java -jar /home/jingjing/Downloads/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=${sample}_aln.sort.bam O=${sample}_aln.rmdup.bam M=${sample}_dup_metrics.txt 
samtools index ${sample}_aln.rmdup.bam

### Separate strand specific alignment_Reverse strand ###
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -b -f 128 -F 16 ${sample}_aln.rmdup.bam > rev1.bam
samtools view -b -f 80 ${sample}_aln.rmdup.bam > rev2.bam

samtools merge -f ${sample}_aln.rmdup.rev.bam rev1.bam rev2.bam
samtools sort ${sample}_aln.rmdup.rev.bam -o ${sample}_aln.rmdup.rev.sort.bam
samtools index ${sample}_aln.rmdup.rev.sort.bam

### Separate strand specific alignment_Forward strand ###

# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand

samtools view -b -f 144 ${sample}_aln.rmdup.bam > fwd1.bam
samtools view -b -f 64 -F 16 ${sample}_aln.rmdup.bam > fwd2.bam

samtools merge -f ${sample}_aln.rmdup.fwd.bam fwd1.bam fwd2.bam
samtools sort ${sample}_aln.rmdup.fwd.bam -o ${sample}_aln.rmdup.fwd.sort.bam
samtools index ${sample}_aln.rmdup.fwd.sort.bam

rm fwd1.bam
rm fwd2.bam
rm rev1.bam
rm rev2.bam
rm ${sample}_aln.rmdup.fwd.bam
rm ${sample}_aln.rmdup.rev.bam
rm ${sample}_aln.map.bam
rm ${sample}_aln.sort.bam
rm ${sample}_aln.sort.bam.bai

bamCoverage --bam ${sample}_aln.rmdup.fwd.sort.bam -o ${sample}_bs5.fwd.bw --normalizeUsing RPKM --binSize 5 -e
bamCoverage --bam ${sample}_aln.rmdup.rev.sort.bam -o ${sample}_bs5.rev.bw --normalizeUsing RPKM --binSize 5 -e


