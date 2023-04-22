#!/bin/bash

sample=pRnhA-1_S4
echo $sample

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



