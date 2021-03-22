########################
#                      #
# Loading the packages #
#                      #
########################

library(ldsep)
library(VariantAnnotation)
library(here)
library(updog)


###################################################################################
#                                                                                 #
# Reading the file and extracting necessary data for genotyping and LD estimation #
#                                                                                 #
###################################################################################

vcfPath <- here('data', 'GWAS_input_data', 'GWAS_geno_input.vcf')
vcfFile <- readVcf(file = vcfPath)
# geno(header(vcfFile)) --> Finding out header names for read depth (DP) and ref allele counts (RO)
sizemat <- geno(vcfFile)$DP
refmat <- geno(vcfFile)$RO


######################################
#                                    #
# Genotyping individuals using updog #
#                                    #
######################################

# This takes quite some time for my data (ca. 6 hours for 35,000 SNPs)
mout <- multidog(refmat = refmat,
                 sizemat = sizemat,
                 ploidy = 4,
                 model = 'norm')

#plot(mout, indices = sample(1:nrow(vcfFile), 3))

# Filter out SNPs that are mostly nonomorphic
msub <- filter_snp(x = mout, pmax(Pr_0, Pr_1, Pr_2, Pr_3, Pr_4) < 0.95)
nrow(msub$snpdf)

# Extracting log-likelihoods of genotypes to be used with 
varnames <- paste0("logL_", 0:4)
varnames
larray <- format_multidog(x = msub, varname = varnames)
class(larray)
dim(larray) # 17178, 135, 5

# Calculation of LD
like_ld <- mldest(geno = larray, K = 4, type = "comp")
