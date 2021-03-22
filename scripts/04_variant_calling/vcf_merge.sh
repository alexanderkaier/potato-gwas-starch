#!/bin/bash

# In- and output data path
IN_PATH=../../data/variant_calls_GWAS

# Compressing and indexing all VCF files
#for vcf in $IN_PATH/*; do
#    bgzip -c $vcf > $vcf.gz
#done

#for vcfgz in $IN_PATH/*.gz; do
#    tabix -p vcf $vcfgz
#done

# tabix -p vcf ../../data/variant_calls_GWAS/freebayes_192_samples_chr01.vcf.gz

# Merging the VCF files from all 12 chromosomes
# Fist, take the header from the first VCF file
grep '^#' $IN_PATH/freebayes_192_samples_chr01.vcf > $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
# Then, concatenate the bodies entries from all files
grep -v '^#' $IN_PATH/freebayes_192_samples_chr01.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr02.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr03.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr04.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr05.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr06.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr07.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr08.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr09.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr10.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr11.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf
grep -v '^#' $IN_PATH/freebayes_192_samples_chr12.vcf >> $IN_PATH/freebayes_192_samples_chr01-chr12.vcf

bgzip -c $IN_PATH/freebayes_192_samples_chr01-chr12.vcf > $IN_PATH/freebayes_192_samples_chr01-chr12.vcf.gz
tabix -p vcf $IN_PATH/freebayes_192_samples_chr01-chr12.vcf.gz