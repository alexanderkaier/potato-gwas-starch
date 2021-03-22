#!/bin/bash

# Defining in- and output files
inAnnotFile=../../data/Reference_genome_and_annotation_file/PGSC_DM_V403_genes_case-corrected.gff
outPredFile=../../data/Post_GWAS/phureja_403_annotoations_refFlat.txt

# Converting the gff3 potato genome annotation file to a refFlat table
gff3ToGenePred $inAnnotFile $outPredFile