#!/bin/bash

# In- and output file paths
IN_PATH=../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth_MissingsSamples-50_SB_BiallSNPs.vcf
OUT_PATH=../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth_MissingsSamples-50_SB_BiallSNPs_AB-0-25.vcf

# Filtering for allele balance
vcffilter -f "! ( AB < 0.25 )" $IN_PATH > $OUT_PATH