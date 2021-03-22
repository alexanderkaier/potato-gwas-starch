#!/bin/bash

IN_PATH=../../data/alignment_data_markeddup/bwa-mem/samtools
# Creating and specifying the directory for the duplicate marking/removal depth statistics
OUT_PATH_MET=../../analysis/alignment_data_markeddup/bwa-mem/samtools

# Number of threads for the system (twice the number of cores, if hyperthreading is supported)
NTHREADS=$(grep -c ^processor /proc/cpuinfo)
# Number of cores of the system
NCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')

for folder in $IN_PATH/*; do
	sample=$(basename $folder)
	baminMarked=${IN_PATH}/${sample}/${sample}_coord_sorted_markeddup.bam
	baminRemoved=${IN_PATH}/${sample}/${sample}_coord_sorted_removeddup.bam
	depthOutMarked=${OUT_PATH_MET}/${sample}/${sample}_markeddup_depth.txt
	depthOutRemoved=${OUT_PATH_MET}/${sample}/${sample}_removeddup_depth.txt

	samtools depth $baminMarked \
		> $depthOutMarked

	samtools depth $baminRemoved \
		> $depthOutRemoved
done