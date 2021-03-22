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
relVcfPathString = "../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30.vcf"
# Turning it into absolute path
vcfPathString = pathlib.Path(relVcfPathString)
absVcfPath = vcfPathString.resolve(strict=True)
absVcfPath = str(absVcfPath)
# Reading it as VCF file
vcfFile = vcf.Reader(filename=absVcfPath)


#####################################################################################
#                                                                                   #
# Reading the mean read depth file for filtering for non-called sites and max depth #
#                                                                                   #
#####################################################################################

# Relative input string for the pandas module
relDepPathString = "../../analysis/variant_calling_GWAS_SNPS/SNPS_135-samples_Biall_SNPs_QUAL-30.ldepth.mean"
# Turning it into absolute path
DepPathString = pathlib.Path(relDepPathString)
absDepPath = DepPathString.resolve(strict=True)
absDepPath = str(absDepPath)
# Reading it as VCF file
meanDepthFrame = pd.read_csv(absDepPath, sep="\t")


#########################################################
#                                                       #
# Creating output file for max. depth filtered variants #
#                                                       #
#########################################################

# Relative input string for the pyVCF module
#relVcfOutPathString = "../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples_QUAL-30_Depth.vcf"
relVcfOutPathString = "../../data/vcf_GWAS_input/freebayes_192_samples_chr01-chr12_BiSNPs_135-samples_QUAL-30_Depth.vcf"
# Turning it into absolute path
vcfOutPathString = pathlib.Path(relVcfOutPathString)
absVcfOutPath = vcfOutPathString.resolve(strict=True)
absVcfOutPath = str(absVcfOutPath)
# Creating the output object using the output string
vcfOutFile = vcf.Writer(open(absVcfOutPath, 'w'), vcfFile)


#####################################################################################
#                                                                                   #
# Actual filtering and writing variants that passed to the  #
#                                                                                   #
#####################################################################################

# Filtering the SNPs according to the formula Depth <= 
'''
for i in range(meanDepthFrame.shape[0]):
    threshold = meanDepthFrame['MEAN_DEPTH'][i] + 10 * meanDepthFrame['MEAN_DEPTH'][i]
    rec = next(vcfFile)
    # Writing selected variant either to stdout or output file, if specified
    if rec.INFO['DP'] >= 20 & rec.INFO['DP'] <= threshold:
        vcfOutFile.write_record(rec)
'''
for i in range(meanDepthFrame.shape[0]):
    threshold = meanDepthFrame['MEAN_DEPTH'][i] + 10 * meanDepthFrame['MEAN_DEPTH'][i]
    rec = next(vcfFile)
    # Writing selected variant to the output file, if specified
    if rec.INFO['DP'] >= 30 & rec.INFO['DP'] <= threshold:
        vcfOutFile.write_record(rec)