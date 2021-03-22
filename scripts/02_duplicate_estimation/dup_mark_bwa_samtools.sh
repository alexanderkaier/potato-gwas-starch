#!/bin/bash

# Program to mark/remove PCR duplicates within the BAM files generated with bwa-mem using samtools markdup


# Specifying the path containing the input data
IN_PATH=../../data/alignment_data_GWAS
# Creaing and specifying the directory for the output alignment files
OUT_PATH_BAM=../../data/alignment_data_markeddup
mkdir -p $OUT_PATH_BAM
# Creating and specifying the directory for the duplicate removal statistics
OUT_PATH_MET=../../analysis/alignment_data_markeddup
mkdir -p $OUT_PATH_MET


# Number of threads for the system (twice the number of cores, if hyperthreading is supported)
NTHREADS=$(grep -c ^processor /proc/cpuinfo)
# Number of cores of the system
NCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')


for folder in $IN_PATH/*; do
	# Extracting the sample name and creating the output folders for each sample
	sample=$(basename $folder)
	mkdir ${OUT_PATH_BAM}/$sample
	mkdir ${OUT_PATH_MET}/$sample


	# Defining the input and output files for the pipeline
	# In- and output for 'samtools sort', name sorted
	# 'samtools sort' automatically appends the '.bam' suffix
	baminSort=${folder}/${sample}_sorted.bam
	bamoutSort=${OUT_PATH_BAM}/${sample}/${sample}_name_sorted.bam
	# Output for 'samtools fixmate'
	bamoutFixmate=${OUT_PATH_BAM}/${sample}/${sample}_name_sorted_fixedmate.bam
	# Output for 'samtools sort', coordinate sorted
	bamoutSortCoord=${OUT_PATH_BAM}/${sample}/${sample}_coord_sorted_fixedmate.bam
	# Output for 'samtools markdup' without removing duplicates
	bamoutMarkeddup=${OUT_PATH_BAM}/${sample}/${sample}_coord_sorted_markeddup.bam
	stderrMarkeddup=${OUT_PATH_MET}/${sample}/${sample}_coord_sorted_markeddup_stderr.stat
	# Output for 'samtools markdup' including removing duplicates
	bamoutRemoveddup=${OUT_PATH_BAM}/${sample}/${sample}_coord_sorted_removeddup.bam
	stderrRemoveddup=${OUT_PATH_MET}/${sample}/${sample}_coord_sorted_removeddup_stderr.stat
	# Output for 'samtools flagstat' without removing duplicates
	flagstatoutMarkeddup=${OUT_PATH_MET}/${sample}/${sample}_flagstat_marked.out
	# Output for 'samtools flagstat' including removing duplicates
	flagstatoutRemoveddup=${OUT_PATH_MET}/${sample}/${sample}_flagstat_removed.out
	# Output for 'samtools stats' without removing duplicates
	samtoolsstatsoutMarkeddup=${OUT_PATH_MET}/${sample}/${sample}_stats_marked.out
	# Output for 'samtools stats' including removing duplicates
	samtoolsstatsoutRemoveddup=${OUT_PATH_MET}/${sample}/${sample}_stats_removed.out
	# Output of computation times for each step
	comptimes=${OUT_PATH_MET}/${sample}/computation_times.txt


	# BAM files need to be name sorted for subsequent mate fixing
	sortn_start=`date +%s`
	samtools sort -@ $NTHREADS -n \
		-o $bamoutSort $baminSort
	sortn_end=`date +%s`
	sortn_time=$((sortn_end-sortn_start))
	echo "samtools sort time in s (name sorted):	$sortn_time" >> $comptimes


	# 'samtools fixmate' fills in mate coordinates. Necessary for 'samtools markdup'
	fixm_start=`date +%s`
	samtools fixmate -@ $NTHREADS -m \
		$bamoutSort $bamoutFixmate
	fixm_end=`date +%s`
	fixm_time=$((fixm_end-fixm_start))
	echo "samtools fixmate time in s:	$fixm_time" >> $comptimes


	# BAM files need to be coordinate sorted for subsequent duplicate marking and removal
	sortc_start=`date +%s`
	samtools sort -@ $NTHREADS \
		-o $bamoutSortCoord $bamoutFixmate
	sortc_end=`date +%s`
	sortc_time=$((sortc_end-sortc_start))
	echo "samtools sort time in s (coordinate sorted):	$sortc_time" >> $comptimes	



	########################################################
	#													   #
	# Creation of BAM files with marked/removed duplicates #
	#													   #
	########################################################

	# Computing the actual PCR duplicates and creating new BAM file with marked diplicates based on the output of 'samtools fixmate'
	md_start=`date +%s`
	samtools markdup -@ $NTHREADS -s \
		$bamoutSortCoord $bamoutMarkeddup \
		2> $stderrMarkeddup
	md_end=`date +%s`
	md_time=$((md_end-md_start))
	echo "samtools markdup time in s (only marking dups):	$md_time" >> $comptimes


	# Computing and removing the actual PCR duplicates and creating new BAM file with marked diplicates based on the output of 'samtools fixmate'
	rm_start=`date +%s`
	samtools markdup -@ $NTHREADS -r -s \
		$bamoutSortCoord $bamoutRemoveddup \
		2> $stderrRemoveddup
	rm_end=`date +%s`
	rm_time=$((rm_end-rm_start))
	echo "samtools markdup time in s (removing dups):	$rm_time" >> $comptimes



	#####################################################
	#													#
	# Indexing BAM files with marked/removed duplicates #
	#												    #
	#####################################################

	# Indexing the BAM file and measuring computation time for marked duplicates
	index_mark_start=`date +%s`
	samtools index $bamoutMarkeddup
	index_mark_end=`date +%s`
	index_mark_time=$((index_mark_end-index_mark_start))
	echo "samtools indexing time in s (only marking dups):	$index_mark_time" >> $comptimes


	# Indexing the BAM file and measuring computation time for removed duplicates
	index_remove_start=`date +%s`
	samtools index $bamoutRemoveddup
	index_remove_end=`date +%s`
	index_remove_time=$((index_remove_end-index_remove_start))
	echo "samtools indexing time in s (removing dups):	$index_remove_time" >> $comptimes



	##################################################################
	#													   			 #
	# Running statistics on BAM files with marked/removed duplicates #
	#													   			 #
	##################################################################

	# Creating stat file using samtools flagstat and measuring computation time for marked duplicates
	flagstat_mark_start=`date +%s`
	samtools flagstat -@ $NTHREADS $bamoutMarkeddup > $flagstatoutMarkeddup
	flagstat_mark_end=`date +%s`
	flagstat_mark_time=$((flagstat_mark_end-flagstat_mark_start))
	echo "samtools flagstat time in s (only marking dups):	$flagstat_mark_time" >> $comptimes


	# Creating stat file using samttols stats and measuring computation time for marked duplicates
	samtools_stats_mark_start=`date +%s`
	samtools stats -@ $NTHREADS $bamoutMarkeddup > $samtoolsstatsoutMarkeddup
	samtools_stats_mark_end=`date +%s`
	samtools_stats_mark_time=$((samtools_stats_mark_end-samtools_stats_mark_start))
	echo "samtools stats time in s (only marking dups):	$samtools_stats_mark_time" >> $comptimes
	

	# Creating stat file using samtools flagstat and measuring computation time for removed duplicates
	flagstat_remove_start=`date +%s`
	samtools flagstat -@ $NTHREADS $bamoutRemoveddup > $flagstatoutRemoveddup
	flagstat_remove_end=`date +%s`
	flagstat_remove_time=$((flagstat_remove_end-flagstat_remove_start))
	echo "samtools flagstat time in s (removing dups):	$flagstat_remove_time" >> $comptimes


	# Creating stat file using samttols stats and measuring computation time for removed duplicates
	samtools_stats_remove_start=`date +%s`
	samtools stats -@ $NTHREADS $bamoutRemoveddup > $samtoolsstatsoutRemoveddup
	samtools_stats_remove_end=`date +%s`
	samtools_stats_remove_time=$((samtools_stats_remove_end-samtools_stats_remove_start))
	echo "samtools stats time in s (removing dups):	$samtools_stats_remove_time" >> $comptimes


	# Recording the total computation time
	echo "total time in s:	$((sortn_time+fixm_time+sortc_time+md_time+rm_time+index_mark_time+index_remove_time+flagstat_mark_time+samtools_stats_mark_time+flagstat_remove_time+samtools_stats_remove_time))" \
		>> $comptimes

done