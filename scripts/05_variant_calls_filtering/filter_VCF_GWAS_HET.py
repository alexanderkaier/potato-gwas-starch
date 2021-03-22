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
relVcfPathString = "../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples_SB.vcf"
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
############################################################

# Relative input string for the pyVCF module
relVcfOutPathString = "../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth_HQGenotypes-MissingSamples_SB_HET.vcf"
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


# Filtering the SNPs according to the minimum number of reference and alternative allele
for rec in vcfFile:
    # Counter for ALT and REF
    REF=0
    ALT=0
    # Checking for the desired qualities in each sample at each site
    for sample in rec.samples:
        if sample.called:
            # Check for existence of '0' and '1' alleles
            if '0' in sample['GT'].split('/'):
                REF+=1
            if '1' in sample['GT'].split('/'):
                ALT+=1
    if REF > 0 and ALT > 0:
        vcfOutFile.write_record(rec)