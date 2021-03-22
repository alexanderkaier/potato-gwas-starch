#!/bin/bash

IN_PATH=../../data/RE_processed
OUT_PATH=../../analysis/raw_seqs
mkdir -p ../../analysis/raw_seqs


# Creating a fastqc quality report for each original bz2 file
#for folder in ${IN_PATH}/*; do
#	sample=$(basename $folder)
#	mkdir -p ${OUT_PATH}/$sample
#	for file in ${folder}/*.bz2; do
#    	fastqc $file -o ${OUT_PATH}/$sample
#    done
#done


# Creating cummulated quality check using multiQC
multiqc $IN_PATH/. -o $OUT_PATH/