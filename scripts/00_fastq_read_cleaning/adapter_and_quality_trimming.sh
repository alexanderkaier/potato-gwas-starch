#!/bin/bash

# Program to filter the raw and compressed fastq.bz2 files by their Q-score and to trim 
# adapters in case of contamination. The arguments in1 and in2 allow analysis of paired reads
# The default kmer handling is kfilter, meaning all observations of a reference kmer will be discarded
# Potential adapter sequences will be trimmed using the built-in adapters.fa file containing all current Illumina adapter sequences

# Path to the reference adapter file and output files
# According to 'https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbduk-guide/',
# this file contains adapters sequences for Illumina TrueSeq and Nextera platforms.
# Whether these adapters are used in SeqNova experiments, was not completely determined.
#ADAPTER_REF=/usr/local/bin/bbmap/resources/adapters.fa
ADAPTER_REF=/home/rna/BioSoftware/bbmap/resources/adapters.fa
PATH_READ_IN=../../data/RE_processed
PATH_READ_OUT=../../data/trimmed_seqs
PATH_STAT_OUT=../../analysis/trimmed_seqs
mkdir -p $PATH_READ_OUT
mkdir -p $PATH_STAT_OUT

for folder in ${PATH_READ_IN}/*; do

        # Creating the directory containing the cleaned FASTQ files and referring to it as output path
        # as well as creating the stat output paths for the statistics and the stdout and stderr
        mkdir ${PATH_READ_OUT}/$(basename $folder)
        mkdir ${PATH_STAT_OUT}/$(basename $folder)
        OUT_PATH=${PATH_READ_OUT}/$(basename $folder)
        STAT_OUT_PATH=${PATH_STAT_OUT}/$(basename $folder)

        R1=$folder/*R1*.bz2
        R1s=$(basename $R1)
        extension_part2=${R1s##*.} # extracting "bz2"
        filename1=${R1s%.*}
        extension_part1=${filename1##*.} # extracting "fastq"
        filename1=${filename1%.*}
        R2=$folder/*R2*.bz2
        R2s=$(basename $R2)
        filename2=${R2s%.*}
        filename2=${filename2%.*}
        extension=$extension_part1.$extension_part2 # creating "fastq.bz2"

        # ktrim (side at which to be trimmed), k (kmer size to be used), hdist (hamming distance/allowed number of mismatches)
        # mink searches for kmers at both ends in the reference and the specified end of the read
        # that are at least of length n in mink=n

        # FIRST, TRIMMING ADAPTERS RIGHT:
        bbduk.sh -Xmx4g \
                in1=${folder}/${R1s} \
                in2=${folder}/${R2s} \
                out1=${OUT_PATH}/${filename1}_clean_right.$extension \
                out2=${OUT_PATH}/${filename2}_clean_right.$extension \
                ref=$ADAPTER_REF \
                ktrim=r \
        	k=23 \
        	mink=11 \
        	hdist=1 \
                stats=${STAT_OUT_PATH}/stats-r.txt \
                > ${STAT_OUT_PATH}/out-r.txt \
                2> ${STAT_OUT_PATH}/err-r.txt


        # SECOND, TRIMMING ADAPTERS LEFT:
        bbduk.sh -Xmx4g \
                in1=${OUT_PATH}/${filename1}_clean_right.$extension \
                in2=${OUT_PATH}/${filename2}_clean_right.$extension \
                out1=${OUT_PATH}/${filename1}_clean_both.$extension \
                out2=${OUT_PATH}/${filename2}_clean_both.$extension \
                ref=$ADAPTER_REF \
                ktrim=l \
        	k=23 \
        	mink=11 \
        	hdist=1 \
                stats=${STAT_OUT_PATH}/stats-b.txt \
                > ${STAT_OUT_PATH}/out-b.txt \
                2> ${STAT_OUT_PATH}/err-b.txt

        # THIRD, TRIMMING BASED ON QUALITY BOTH SIDES AT THE SAME TIME:
        bbduk.sh -Xmx4g \
                in1=${OUT_PATH}/${filename1}_clean_both.$extension \
                in2=${OUT_PATH}/${filename2}_clean_both.$extension \
                out1=${OUT_PATH}/${filename1}_clean30.$extension \
                out2=${OUT_PATH}/${filename2}_clean30.$extension \
                qtrim=rl \
                trimq=30 \
                minlen=35 \
                stats=${STAT_OUT_PATH}/stats-c.txt \
                > ${STAT_OUT_PATH}/out-c.txt \
                2> ${STAT_OUT_PATH}/err-c.txt

        # Deleting intermediate FASTQ files, if desired
        #rm ${OUT_PATH}/*right.$extension ${OUT_PATH}/*both.$extension
done