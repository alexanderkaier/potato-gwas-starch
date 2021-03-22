#!/bin/bash

IN_PATH=../../data/alignment_data
# Creating and specifying the directory for the duplicate marking/removal depth statistics
OUT_PATH=../../analysis/alignment_data

# Number of threads for the system (twice the number of cores, if hyperthreading is supported)
NTHREADS=$(grep -c ^processor /proc/cpuinfo)
# Number of cores of the system
NCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')

# Looping over all samples within both aligner folders and calculating the depth at every covered position
for alignerFolder in $IN_PATH/*; do
	for sampleFolder in $alignerFolder/*; do
		aligner=$(basename $alignerFolder)
		sample=$(basename $sampleFolder)
		bamIn=${IN_PATH}/${aligner}/${sample}/${sample}_sorted.bam
		depthOut=${OUT_PATH}/${aligner}/${sample}/${sample}_depth.txt
		
		samtools depth $bamIn > $depthOut
	done
done