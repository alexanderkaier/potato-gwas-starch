#!/bin/bash

# In- and output file paths
IN_PATH=../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth.vcf
OUT_PATH=../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth_HQ-Genotypes.vcf

# Filtering for allele bias of alternative alleles
vcffilter -g "DP > 60 & GQ > 20" $IN_PATH > $OUT_PATH