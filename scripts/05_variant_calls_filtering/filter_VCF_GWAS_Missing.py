#!/bin/python

# importing modules
import pathlib
import vcf
import pandas as pd
import numpy as np

##############################################################################
#                                                                            #
# Reading the VCF file for filtering out all sites with no genotyped samples #
#                                                                            #
##############################################################################

# Relative input string for the pyVCF module
relVcfPathString = "../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth.vcf"
# Turning it into absolute path
vcfPathString = pathlib.Path(relVcfPathString)
absVcfPath = vcfPathString.resolve(strict=True)
absVcfPath = str(absVcfPath)
# Reading it as VCF file
vcfFile = vcf.Reader(filename=absVcfPath)

#############################################################
#                                                           #
# Creating output file for missing sample filtered variants #
#                                                           #
#############################################################

# Relative input string for the pyVCF module
#relVcfOutPathString = "../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth_MissingsSamples-75.vcf"
relVcfOutPathString = "../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples.vcf"
# Turning it into absolute path
vcfOutPathString = pathlib.Path(relVcfOutPathString)
absVcfOutPath = vcfOutPathString.resolve(strict=True)
absVcfOutPath = str(absVcfOutPath)
# Creating the output object using the output string
vcfOutFile = vcf.Writer(open(absVcfOutPath, 'w'), vcfFile)

#############################################################################################
#                                                                                           #
# Actual filtering and writing variants that passed the missing sample percentage threshold #
#                                                                                           #
#############################################################################################

# Filtering the SNPs according to missing data of 10% (i.e. tolerating up to 90% missing individuals at any site)
for rec in vcfFile:
    # Checking for the desired qualities in each sample at each site
    for sample in rec.samples:
        if sample.called:
            if sample['DP'] < 61 or sample['GQ'] < 20:
                sample.called = False
            elif '1' in sample['GT'].split('/') and sample['AO']/sample['DP'] < 0.25:
                    sample.called = False

    # Filtering for call rate
    if rec.call_rate >= 0.90:
        vcfOutFile.write_record(rec)