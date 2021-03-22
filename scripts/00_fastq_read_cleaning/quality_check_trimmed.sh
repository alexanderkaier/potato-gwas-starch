#!/bin/bash

IN_PATH=../../data/trimmed_seqs
OUT_PATH=../../analysis/trimmed_seqs/
mkdir -p ../../analysis/trimmed_seqs

# Creating a fastqc quality report for each trimmed bz2 file
for folder in ${IN_PATH}/*; do
	SAMPLE=$(basename $folder)
	for file in $folder/*clean30.fastq.bz2; do
    	fastqc $file -o ${OUT_PATH}$SAMPLE
	done
done

# Creating cummulated quality check using multiQC
multiqc $OUT_PATH -o $OUT_PATH