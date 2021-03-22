#!/bin/bash

###################
#                 #
# In- and outputs #
#                 #
###################

# In- and output for subsetting the SNPs variant call 
IN_PATH_SUBSET=../../data/variant_calls_GWAS_filtered

# In- and output for the statistics
:'
IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_SNPs_subset100000.vcf
OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS_subset100000
mkdir -p $OUT_STAT_PATH
OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_subset100000

# Subsetting ca. 100,000 random SNPs from the VCF file containing only SNPs
bcftools view $IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_SNPs.vcf | vcfrandomsample --rate 0.054 \
    > $IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_SNPs_subset100000.vcf


# Extracting information from the subsetted SNP VCF file
#vcftools --vcf $IN_VCF_FILE --freq2 --out $OUT_STAT_NAME --max-alleles 4 --> vcftools does not work on ploidy-dependent stuff :(
vcftools --vcf $IN_VCF_FILE --site-mean-depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME
#vcftools --vcf $IN_VCF_FILE --missing-indv --out $OUT_STAT_NAME --> vcftools does not work on ploidy-dependent stuff :(
vcftools --vcf $IN_VCF_FILE --missing-site --out $OUT_STAT_NAME
#vcftools --vcf $IN_VCF_FILE --het --out $OUT_STAT_NAME --> vcftools does not work on ploidy-dependent stuff :(
'

IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_SNPs_135-samples_subset100000.vcf
OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS_subset100000
mkdir -p $OUT_STAT_PATH
OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_135-samples_subset100000


# Extracting information from the subsetted SNP VCF file
#vcftools --vcf $IN_VCF_FILE --freq2 --out $OUT_STAT_NAME --max-alleles 4 --> vcftools does not work on ploidy-dependent stuff :(
vcftools --vcf $IN_VCF_FILE --site-mean-depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME
#vcftools --vcf $IN_VCF_FILE --missing-indv --out $OUT_STAT_NAME --> vcftools does not work on ploidy-dependent stuff :(
#vcftools --vcf $IN_VCF_FILE --missing-site --out $OUT_STAT_NAME --> vcftools does not work on ploidy-dependent stuff :(
#vcftools --vcf $IN_VCF_FILE --het --out $OUT_STAT_NAME --> vcftools does not work on ploidy-dependent stuff :(
