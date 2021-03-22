#!/bin/bash

# In- and output file paths
IN_PATH=../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples.vcf
OUT_PATH=../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples_SB.vcf

# Filtering for allele bias of alternative alleles
vcffilter -f "SAF > 0 & SAR > 0" $IN_PATH > $OUT_PATH