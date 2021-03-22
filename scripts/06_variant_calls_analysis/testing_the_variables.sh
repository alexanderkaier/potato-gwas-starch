#!/bin/bash

IN_PATH=../../data/variant_calls

# Creating all output paths for filter strategies
PATH_OUT_VCF=../../data/variant_calls_filtered

for vcfFile in $IN_PATH/*.vcf; do
    # Extracting the file name for subsequent appending
    file=$(basename $vcfFile)
    extension=${file##*.} # extracting "vcf"
    filename=${file%.*} # Extracting the name w\o extension
    echo $file
    echo $filename
    echo $extension

    echo "${filename}_filtered_.$extension"
done
