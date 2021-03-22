#!/bin/bash

# Program to mark/remove PCR duplicates within the BAM files generated with bwa-mem using samtools markdup


# Specifying the path containing the input data (alignments without quality-trimming reads)
IN_PATH=../../data/alignment_data_duplicates
# Creaing and specifying the directory for the output alignment files
OUT_PATH_BAM=../../data/alignment_data_markeddup
mkdir -p $OUT_PATH_BAM
# Creating and specifying the directory for the duplicate removal statistics
OUT_PATH_MET=../../analysis/alignment_data_markeddup
mkdir -p  $OUT_PATH_MET 


# Number of threads for the system (twice the number of cores, if hyperthreading is supported)
NTHREADS=$(grep -c ^processor /proc/cpuinfo)
# Number of cores of the system
NCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')

# Initializing counter for stdout progress messages
count=1
# Calculating number of total input samples
TOTAL=$(ls -l $IN_PATH | grep -c ^d)

for folder in $IN_PATH/*; do
	# Extracting the sample name and creating the output folders for each sample
	sample=$(basename $folder)
	mkdir -p ${OUT_PATH_BAM}/$sample
	mkdir -p ${OUT_PATH_MET}/$sample


	# Defining the input and output files for the pipeline
	# In- and output for 'samtools sort', name sorted
	# 'samtools sort' automatically appends the '.bam' suffix
	BAM_IN_SORT_NAME=${folder}/${sample}_sorted.bam
	BAM_OUT_SORT_NAME=${OUT_PATH_BAM}/${sample}/${sample}_name_sorted.bam
	# Output for 'samtools fixmate'
	BAM_OUT_FIXMATE=${OUT_PATH_BAM}/${sample}/${sample}_name_sorted_fixedmate.bam
	# Output for 'samtools sort', coordinate sorted
	BAM_OUT_SORT_COORD=${OUT_PATH_BAM}/${sample}/${sample}_coord_sorted_fixedmate.bam
	# Output for 'samtools markdup' without removing duplicates
	BAM_OUT_MARKEDDUP=${OUT_PATH_BAM}/${sample}/${sample}_coord_sorted_markeddup.bam
	STDERR_MARKEDDUP=${OUT_PATH_MET}/${sample}/${sample}_coord_sorted_markeddup_stderr.stat
	# Output for 'samtools markdup' including removing duplicates
	BAM_OUT_REMOVEDDUP=${OUT_PATH_BAM}/${sample}/${sample}_coord_sorted_removeddup.bam
	STDERR_REMOVEDDUP=${OUT_PATH_MET}/${sample}/${sample}_coord_sorted_removeddup_stderr.stat
	# Output for 'samtools flagstat' without removing duplicates
	FLAGSTAT_OUT_MARKEDDUP=${OUT_PATH_MET}/${sample}/${sample}_flagstat_marked.out
	# Output for 'samtools flagstat' including removing duplicates
	FLAGSTAT_OUT_REMOVEDDUP=${OUT_PATH_MET}/${sample}/${sample}_flagstat_removed.out
	# Output for 'samtools stats' without removing duplicates
	STATS_OUT_MARKEDDUP=${OUT_PATH_MET}/${sample}/${sample}_stats_marked.out
	# Output for 'samtools stats' including removing duplicates
	STATS_OUT_REMOVEDDUP=${OUT_PATH_MET}/${sample}/${sample}_stats_removed.out
	# Output of computation times for each step
	COMP=${OUT_PATH_MET}/${sample}/computation_times.txt


	# BAM files need to be name sorted for subsequent mate fixing
	sortn_start=`date +%s`
	samtools sort -@ $NTHREADS -n \
		-o $BAM_OUT_SORT_NAME $BAM_IN_SORT_NAME
	sortn_end=`date +%s`
	sortn_time=$((sortn_end-sortn_start))
	echo "samtools sort time in s (name sorted):	$sortn_time" >> $COMP
	echo "Name sorting finished for sample $sample ($count of $TOTAL)"


	# 'samtools fixmate' fills in mate coordinates. Necessary for 'samtools markdup'
	fixm_start=`date +%s`
	samtools fixmate -@ $NTHREADS -m \
		$BAM_OUT_SORT_NAME $BAM_OUT_FIXMATE
	fixm_end=`date +%s`
	fixm_time=$((fixm_end-fixm_start))
	echo "samtools fixmate time in s:	$fixm_time" >> $COMP
	echo "Mate fixing finished for sample $sample ($count of $TOTAL)"


	# BAM files need to be coordinate sorted for subsequent duplicate marking and removal
	sortc_start=`date +%s`
	samtools sort -@ $NTHREADS \
		-o $BAM_OUT_SORT_COORD $BAM_OUT_FIXMATE
	sortc_end=`date +%s`
	sortc_time=$((sortc_end-sortc_start))
	echo "samtools sort time in s (coordinate sorted):	$sortc_time" >> $COMP
	echo "Coordinate sorting finished for sample $sample ($count of $TOTAL)"	



	########################################################
	#													   #
	# Creation of BAM files with marked/removed duplicates #
	#													   #
	########################################################

	# Computing the actual PCR duplicates and creating new BAM file with marked duplicates based on the output of 'samtools fixmate'
	md_start=`date +%s`
	samtools markdup -@ $NTHREADS -s \
		$BAM_OUT_SORT_COORD $BAM_OUT_MARKEDDUP \
		2> $STDERR_MARKEDDUP
	md_end=`date +%s`
	md_time=$((md_end-md_start))
	echo "samtools markdup time in s (only marking dups):	$md_time" >> $COMP
	echo "Duplicate marking finished for sample $sample ($count of $TOTAL)"


	# Computing and removing the actual PCR duplicates and creating new BAM file with marked diplicates based on the output of 'samtools fixmate'
	rm_start=`date +%s`
	samtools markdup -@ $NTHREADS -r -s \
		$BAM_OUT_SORT_COORD $BAM_OUT_REMOVEDDUP \
		2> $STDERR_REMOVEDDUP
	rm_end=`date +%s`
	rm_time=$((rm_end-rm_start))
	echo "samtools markdup time in s (removing dups):	$rm_time" >> $COMP
	echo "Duplicate removal finished for sample $sample ($count of $TOTAL)"



	#####################################################
	#													#
	# Indexing BAM files with marked/removed duplicates #
	#												    #
	#####################################################

	# Indexing the BAM file and measuring computation time for marked duplicates
	index_mark_start=`date +%s`
	samtools index $BAM_OUT_MARKEDDUP
	index_mark_end=`date +%s`
	index_mark_time=$((index_mark_end-index_mark_start))
	echo "samtools indexing time in s (only marking dups):	$index_mark_time" >> $COMP
	echo "Indexing of BAM with marked duplicates finished for sample $sample ($count of $TOTAL)"


	# Indexing the BAM file and measuring computation time for removed duplicates
	index_remove_start=`date +%s`
	samtools index $BAM_OUT_REMOVEDDUP
	index_remove_end=`date +%s`
	index_remove_time=$((index_remove_end-index_remove_start))
	echo "samtools indexing time in s (removing dups):	$index_remove_time" >> $COMP
	echo "Indexing of BAM with removed duplicates finished for sample $sample ($count of $TOTAL)"



	##################################################################
	#													   			 #
	# Running statistics on BAM files with marked/removed duplicates #
	#													   			 #
	##################################################################

	# Creating stat file using samtools flagstat and measuring computation time for marked duplicates
	flagstat_mark_start=`date +%s`
	samtools flagstat -@ $NTHREADS $BAM_OUT_MARKEDDUP > $FLAGSTAT_OUT_MARKEDDUP
	flagstat_mark_end=`date +%s`
	flagstat_mark_time=$((flagstat_mark_end-flagstat_mark_start))
	echo "samtools flagstat time in s (only marking dups):	$flagstat_mark_time" >> $COMP
	echo "Summary creation with 'samtools flagstat' of BAM with marked duplicates finished for sample $sample ($count of $TOTAL)"


	# Creating stat file using samttols stats and measuring computation time for marked duplicates
	samtools_stats_mark_start=`date +%s`
	samtools stats -@ $NTHREADS $BAM_OUT_MARKEDDUP > $STATS_OUT_MARKEDDUP
	samtools_stats_mark_end=`date +%s`
	samtools_stats_mark_time=$((samtools_stats_mark_end-samtools_stats_mark_start))
	echo "samtools stats time in s (only marking dups):	$samtools_stats_mark_time" >> $COMP
	echo "Summary creation with 'samtools stats' of BAM with marked duplicates finished for sample $sample ($count of $TOTAL)"
	

	# Creating stat file using samtools flagstat and measuring computation time for removed duplicates
	flagstat_remove_start=`date +%s`
	samtools flagstat -@ $NTHREADS $BAM_OUT_REMOVEDDUP > $FLAGSTAT_OUT_REMOVEDDUP
	flagstat_remove_end=`date +%s`
	flagstat_remove_time=$((flagstat_remove_end-flagstat_remove_start))
	echo "samtools flagstat time in s (removing dups):	$flagstat_remove_time" >> $COMP
	echo "Summary creation with 'samtools flagstat' of BAM with removed duplicates finished for sample $sample ($count of $TOTAL)"


	# Creating stat file using samttols stats and measuring computation time for removed duplicates
	samtools_stats_remove_start=`date +%s`
	samtools stats -@ $NTHREADS $BAM_OUT_REMOVEDDUP > $STATS_OUT_REMOVEDDUP
	samtools_stats_remove_end=`date +%s`
	samtools_stats_remove_time=$((samtools_stats_remove_end-samtools_stats_remove_start))
	echo "samtools stats time in s (removing dups):	$samtools_stats_remove_time" >> $COMP
	echo "Summary creation with 'samtools stats' of BAM with removed duplicates finished for sample $sample ($count of $TOTAL)"


	# Recording the total computation time
	echo "total time in s:	$((sortn_time+fixm_time+sortc_time+md_time+rm_time+index_mark_time+index_remove_time+flagstat_mark_time+samtools_stats_mark_time+flagstat_remove_time+samtools_stats_remove_time))" \
		>> $COMP
	echo "Sample $sample done ($count of $TOTAL)!"

	# Incrementing the counter of samples processed
	(( count++ ))

done