#!/bin/bash


# In- and output for vcffilter
IN_VCF_FILE=../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples.vcf
OUT_VCF_FILE=../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30.vcf

# Filtering all sites with QUAL < 30
vcffilter -f " ! ( QUAL < 30 ) " \
    $IN_VCF_FILE > $OUT_VCF_FILE