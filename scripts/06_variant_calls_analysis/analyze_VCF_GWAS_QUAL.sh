#!/bin/bash

# Checking for different metrics for _next filtering_step
IN_PATH_SUBSET=../../data/vcf_GWAS_input
IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples_SB_HET.vcf
OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS
mkdir -p $OUT_STAT_PATH
OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_GWAS


# Extracting information from the subsetted SNP VCF file
vcftools --vcf $IN_VCF_FILE --site-mean-depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME
