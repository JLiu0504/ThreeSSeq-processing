#!/bin/bash

for read1 in *_L001_R1_001.fastq 

do
	echo $read1
        output=${read1/L001_R1_001.fastq/R1.trim.fastq}
	echo $output
        java -jar /home/jingjing/Downloads/Trimmomatic-0.36/trimmomatic-0.36.jar SE $read1 $output HEADCROP:3
	
done


