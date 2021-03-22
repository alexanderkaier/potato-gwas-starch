#!/bin/python3

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

###############################################
#                                             #
# Creating output file for annotated variants #
#                                             #
###############################################

# Relative input string for the pyVCF module
relVcfOutPathString = "../../data/GWAS_input_data/GWAS_geno_input.vcf"
# Turning it into absolute path
vcfOutPathString = pathlib.Path(relVcfOutPathString)
absVcfOutPath = vcfOutPathString.resolve(strict=True)
absVcfOutPath = str(absVcfOutPath)
# Creating the output object using the output string
vcfOutFile = vcf.Writer(open(absVcfOutPath, 'w'), vcfFile)

# Creating IDs for each variant comprising of CHROM name and POS on that choromosome
for rec in vcfFile:
    rec.ID = str(rec.CHROM) + '_' + str(rec.POS)
    vcfOutFile.write_record(rec)
