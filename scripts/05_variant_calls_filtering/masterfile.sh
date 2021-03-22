#!/bin/bash

################################################################
#                                                              #
# Master script for applying all selected filters sequencially #
#                                                              #
################################################################

: '
We start with the VCF file only containing SNPs and 
phenotyped samples
'


# Only biallelic SNPs
touch ../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples.vcf
python3 filter_VCF_GWAS_biallelic_SNPs.py

# QUAL >= 30
sh filter_VCF_GWAS_QUAL.sh
# Create statistics for average depth for subsequent filtering of max depth
IN_PATH_SUBSET=../../data/vcf_GWAS_input
IN_VCF_FILE=$IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30.vcf
OUT_STAT_PATH=../../analysis/variant_calling_GWAS_SNPS
mkdir -p $OUT_STAT_PATH
OUT_STAT_NAME=$OUT_STAT_PATH/SNPS_135-samples_Biall_SNPs_QUAL-30
# Extracting information from the subsetted SNP VCF file
vcftools --vcf $IN_VCF_FILE --site-mean-depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --depth --out $OUT_STAT_NAME
vcftools --vcf $IN_VCF_FILE --site-quality --out $OUT_STAT_NAME


# Depth filter
touch ../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth.vcf
python3 filter_VCF_GWAS_depth_per_site.py


# Blanking "bad" genotypes and filtering variants with too few genotyped samples 
touch ../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples.vcf
python3 filter_VCF_GWAS_Missing.py


# Min. 1 read per strand per variant
sh filter_VCF_GWAS_StrandBias.sh


# MAF 0.5 %
#sh filter_VCF_GWAS_MAF.sh


# Min. 1 ALT and 1 REF allele in the population at this site
touch ../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples_SB_HET.vcf
python3 filter_VCF_GWAS_HET.py