#!/bin/bash

# This script creates and index of the reference file for alignments using the bwa suit. This has to be done only once.
# WARNING: This file has yet to be modified to store the the produced index by 'bwa index' in another folder

# Specifying the directory of the reference genome
REF_GENOME=../../data/Reference_genome_and_annotation_file/potato_dm_v404_all_pm_un.fasta
# Creating the index
bwa index $REF_GENOME

