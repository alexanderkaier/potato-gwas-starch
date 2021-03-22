#!/bin/bash

###############################################################################################################################################
#																																  		      #
# Program for extracting compressed FASTQ files, performing alignment using bwa-mem on adapter-trimmed and quality-filtered reads,		      #
# creating a statistics output for every alignment, converting alignment directly to BAM without generation of a SAM file,		 		      #
# as well as sorting and indexing all created BAM files.																		  		      #
# Additionally, the time for each step as well as the total time is recorded for every sample and a short message is printed to the terminal  #
#																												 				  		      #
###############################################################################################################################################

# Specifying the directory of the reference genome
REF_GENOME=../../data/Reference_genome_and_annotation_file/potato_dm_v404_all_pm_un.fasta
# Specifying the path containing the input data.
IN_PATH=../../data/trimmed_seqs
# Creaing and specifying the directory for the output alignment files
OUT_PATH_BAM=../../data/alignment_data_GWAS
mkdir -p $OUT_PATH_BAM
# Creating and specifying the directory for the alignment statistics
OUT_PATH_STAT=../../analysis/alignment_data_GWAS
mkdir -p $OUT_PATH_STAT


# Number of threads for the system (twice the number of cores, if hyperthreading is supported)
NTHREADS=$(grep -c ^processor /proc/cpuinfo)
# Number of cores of the system
NCORES=$(grep ^cpu\\scores /proc/cpuinfo | uniq |  awk '{print $4}')


# Initializing counter for stdout progress messages
count=1
# Calculating number of total input samples
TOTAL=$(ls -l $IN_PATH | grep -c ^d)

for folder in $IN_PATH/*; do
	# Specifying the ouput name and creating the output folders for each sample
	sample=$(basename $folder)
	mkdir ${OUT_PATH_BAM}/$sample
	mkdir ${OUT_PATH_STAT}/$sample
	# Defining variables
	INIT_BAM=${OUT_PATH_BAM}/${sample}/${sample}.bam
	SORT_BAM=${OUT_PATH_BAM}/${sample}/${sample}_sorted.bam
	COMP=${OUT_PATH_STAT}/${sample}/computation_times.txt


	# Extracting the compressed FASTQ files of the adapter-trimmed sequences for alignment and measuring computation time
	unzip_start=`date +%s`
	for comp in ${folder}/*clean30.fastq.bz2; do
		bzip2 -dk $comp
	done
	unzip_end=`date +%s`
	unzip_time=$((unzip_end-unzip_start))
	echo "FASTQ extraction time in s:	$unzip_time" >> $COMP
	echo "Extracted both FASTQ files of sample $sample ($count of $TOTAL)"


	# Specifying the input files for alignment
	fastq1=${folder}/*R1*.fastq
	fastq2=${folder}/*R2*.fastq


	# Conducting alignment and measuring computation time
	align_start=`date +%s`
	bwa mem  -t $NTHREADS $REF_GENOME $fastq1 $fastq2 \
		|samtools view -@ $NTHREADS -ub - > $INIT_BAM
	align_end=`date +%s`
	align_time=$((align_end-align_start))
	echo "bwa mem time in s:	$align_time" >> $COMP
	echo "Alignment and compression to BAM file done for sample ($count of $TOTAL)"


	# Creating stat file using samtools flagstat and measuring computation time
	flagstat_start=`date +%s`
	samtools flagstat -@ $NTHREADS $INIT_BAM > ${OUT_PATH_STAT}/${sample}/flagstat_stats.out
	flagstat_end=`date +%s`
	flagstat_time=$((flagstat_end-flagstat_start))
	echo "samtools flagstat time in s:	$flagstat_time" >> $COMP
	echo "samtools flagstat finished for sample $sample ($count of $TOTAL)"


	# Sorting the output BAM file after coordinate and measuring computation time
	sort_start=`date +%s`
	samtools sort -@ $NTHREADS -o $SORT_BAM $INIT_BAM
	sort_end=`date +%s`
	sort_time=$((sort_end-sort_start))
	echo "samtools sorting time in s:	$sort_time" >> $COMP
	echo "samtools sorting finished for sample $sample ($count of $TOTAL)"


	# Creating stat file using samttols stats and measuring computation time
	samtools_stats_start=`date +%s`
	samtools stats -@ $NTHREADS $SORT_BAM > ${OUT_PATH_STAT}/${sample}/stats_stats.out
	samtools_stats_end=`date +%s`
	samtools_stats_time=$((samtools_stats_end-samtools_stats_start))
	echo "samtools stats time in s:	$samtools_stats_time" >> $COMP
	echo "samtools stats finished for sample $sample ($count of $TOTAL)"


	# Indexing the sorted BAM file and measuring computation time
	index_start=`date +%s`
	samtools index $SORT_BAM
	index_end=`date +%s`
	index_time=$((index_end-index_start))
	echo "samtools indexing time in s:	$index_time" >> $COMP
	echo "samtools indexing finished for $sample ($count of $TOTAL)"


	# Running stats on the indexed file and measuring computation time
	idxstats_start=`date +%s`
	samtools idxstats $SORT_BAM > ${OUT_PATH_STAT}/${sample}/idxstat_stats.out
	idxstats_end=`date +%s`
	idxstats_time=$((idxstats_end-idxstats_start))
	echo "samtools idxstats time in s:	$idxstats_time" >> $COMP
	echo "samtools idxstats finished for $sample ($count of $TOTAL)"


	# Recording the total computation time
	echo "total time in s:	$((unzip_time+align_time+flagstat_time+sort_time+samtools_stats_time+index_time+idxstats_time))" \
		>> $COMP
	echo "Sample $sample done ($count of $TOTAL)!"

	# Incrementing the counter of samples processed
	(( count++ ))

	# Deleting uncompressed FASTQ files from input directory
	for fastq in ${folder}/*.fastq; do
		rm $fastq
	done
done
