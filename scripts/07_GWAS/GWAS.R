###############################
#                             #
# Load the necessary packages #
#                             #
###############################

library(here)
library(GWASpoly)
library(gemma2)
library(tidyverse)
library(heritability)
library(cowplot)
library(qqman)
library(ldsep)


#########################################
#                                       #
# Reading the geno- and phenotype files #
#                                       #
#########################################

# Read the genotype and phenotype csv files
genoFile <- here('data', 'GWAS_input_data', 'GWAS_geno_input.csv')
phenoFile <- here('data', 'GWAS_input_data', 'GWAS_pheno_input.csv')
phenoData <- read_csv(phenoFile)
ggplot(phenoData, aes(x=Starch_ratio)) + 
  geom_histogram()
mean(phenoData$Starch_ratio); mean(phenoData$Starch_ratio) + 2 * sd(phenoData$Starch_ratio); mean(phenoData$Starch_ratio) - 2 * sd(phenoData$Starch_ratio)


data <- read.GWASpoly(ploidy = 4, pheno.file = phenoFile,
                      geno.file = genoFile, format = "ACGT", n.traits = 1, delim = ',')


###############################################
#                                             #
# Creating kinship matrix and performing GWAS #
#                                             #
###############################################

# Set the kinship matrix for polygenic effects
data <- set.K(data)

# Further incorporation of population structure
params <- set.params(MAF = 0.01, P3D = F)

#Performing the GWAS multiple times
dataAdd <- GWASpoly(data = data, models = c("additive"),
                  traits = "Starch_ratio", params = params)
dataOne <- GWASpoly(data = data, models = c("1-dom"),
                  traits = "Starch_ratio", params = params)
dataTwo <- GWASpoly(data = data, models = c("2-dom"),
                  traits = "Starch_ratio", params = params)
dataDipGen <- GWASpoly(data = data, models = c("diplo-general"),
                       traits = "Starch_ratio", params = params)
dataDipAdd <- GWASpoly(data = data, models = c("diplo-additive"),
                       traits = "Starch_ratio", params = params)
dataGen <- GWASpoly(data = data, models = c("general"),
                       traits = "Starch_ratio", params = params)


####################################################################################
#                                                                                  #
# Heritability estimation, lambda estimation and GWAS sig. threshold determination #
#                                                                                  #
####################################################################################

# Estimating narrow-sense heritability for all models
# Additive model
h2AddEst <- marker_h2(data.vector = dataAdd@pheno$Starch_ratio, geno.vector = dataAdd@pheno$Individual,
                  K = dataAdd@K, max.iter = 1000); h2Add <- h2AddEst$h2; h2AddInt <- h2AddEst$conf.int1
# Simplex dominant
h2OneEst <- marker_h2(data.vector = dataOne@pheno$Starch_ratio, geno.vector = dataOne@pheno$Individual,
                   K = dataOne@K, max.iter = 1000); h2One <- h2OneEst$h2; h2OneInt <- h2OneEst$conf.int1
# Duplex dominant
h2TwoEst <- marker_h2(data.vector = dataTwo@pheno$Starch_ratio, geno.vector = dataTwo@pheno$Individual,
                   K = dataTwo@K, max.iter = 1000); h2Two <- h2TwoEst$h2; h2TwoInt <- h2TwoEst$conf.int1
# Diploid General
h2DipGenEst <- marker_h2(data.vector = dataDipGen@pheno$Starch_ratio, geno.vector = dataDipGen@pheno$Individual,
                      K = dataDipGen@K, max.iter = 1000); h2DipGen <- h2DipGenEst$h2; h2DipGenInt <- h2DipGenEst$conf.int1
# Diploid Additive
h2DipAddEst <- marker_h2(data.vector = dataDipAdd@pheno$Starch_ratio, geno.vector = dataDipAdd@pheno$Individual,
                      K = dataDipAdd@K, max.iter = 1000); h2DipAdd <- h2DipAddEst$h2; h2DipAddInt <- h2DipAddEst$conf.int1
# General
h2GenEst <- marker_h2(data.vector = dataGen@pheno$Starch_ratio, geno.vector = dataGen@pheno$Individual,
                      K = dataGen@K, max.iter = 1000); h2Gen <- h2GenEst$h2; h2GenInt <- h2GenEst$conf.int1

# Setting significance thresholds
dataAdd <- set.threshold(data = dataAdd, method = "M.eff", level = 0.05)
dataAdd <- set.threshold(data = dataAdd, method = "Bonferroni", level = 0.05)
dataAdd <- set.threshold(data = dataAdd, method = "FDR", level = 0.05)
dataAdd <- set.threshold(data = dataAdd, method = "permute", level = 0.05, n.permute = 100)

dataOne <- set.threshold(data = dataOne, method = "M.eff", level = 0.05)
dataOne <- set.threshold(data = dataOne, method = "Bonferroni", level = 0.05)
dataOne <- set.threshold(data = dataOne, method = "FDR", level = 0.05)
dataOne <- set.threshold(data = dataOne, method = "permute", level = 0.05, n.permute = 100)

dataTwo <- set.threshold(data = dataTwo, method = "M.eff", level = 0.05)
dataTwo <- set.threshold(data = dataTwo, method = "Bonferroni", level = 0.05)
dataTwo <- set.threshold(data = dataTwo, method = "FDR", level = 0.05)
dataTwo <- set.threshold(data = dataTwo, method = "permute", level = 0.05, n.permute = 100)

dataDipGen <- set.threshold(data = dataDipGen, method = "M.eff", level = 0.05)
dataDipGen <- set.threshold(data = dataDipGen, method = "Bonferroni", level = 0.05)
dataDipGen <- set.threshold(data = dataDipGen, method = "FDR", level = 0.05)
dataDipGen <- set.threshold(data = dataDipGen, method = "permute", level = 0.05, n.permute = 100)

dataDipAdd <- set.threshold(data = dataDipAdd, method = "M.eff", level = 0.05)
dataDipAdd <- set.threshold(data = dataDipAdd, method = "Bonferroni", level = 0.05)
dataDipAdd <- set.threshold(data = dataDipAdd, method = "FDR", level = 0.05)
dataDipAdd <- set.threshold(data = dataDipAdd, method = "permute", level = 0.05, n.permute = 100)

dataGen <- set.threshold(data = dataGen, method = "M.eff", level = 0.05)
dataGen <- set.threshold(data = dataGen, method = "Bonferroni", level = 0.05)
dataGen <- set.threshold(data = dataGen, method = "FDR", level = 0.05)
dataGen <- set.threshold(data = dataGen, method = "permute", level = 0.05, n.permute = 100)

qqAdd <- qq.plot(dataAdd) + ggtitle("Additive K model") + 
  xlim(0, 5) + 
  ylim(0, 6.5) + 
  geom_hline(aes(yintercept=4.91), color = "#1B9E77")
qqOne <- qq.plot(dataOne) + ggtitle("Simplex K model")

qqTwo <- qq.plot(dataTwo) + ggtitle("Duplex K model")

qqDipGen <- qq.plot(dataDipGen) + ggtitle(" Diploid general K model")

qqDipAdd <- qq.plot(dataDipAdd) + ggtitle("Diploid additive K model")

qqGen <- qq.plot(dataGen) + ggtitle("General K model")











# write results to the output files
# Additive model
outFileScoresAdd <- here('data', 'GWAS_results', 'scoresAdd.csv')
write.GWASpoly(data = dataAdd, trait = 'Starch_ratio', filename = outFileScoresAdd, what = 'scores', delim = ',')
outFileEffectsAdd <- here('data', 'GWAS_results', 'effectsAdd.csv')
write.GWASpoly(data = dataAdd, trait = 'Starch_ratio', filename = outFileEffectsAdd, what = 'effects', delim = ',')
# Simplex dominant model
outFileScoresOne <- here('data', 'GWAS_results', 'scoresOne.csv')
write.GWASpoly(data = dataOne, trait = 'Starch_ratio', filename = outFileScoresOne, what = 'scores', delim = ',')
outFileEffectsOne <- here('data', 'GWAS_results', 'effectsOne.csv')
write.GWASpoly(data = dataOne, trait = 'Starch_ratio', filename = outFileEffectsOne, what = 'effects', delim = ',')
# Duplex dominant model
outFileScoresTwo <- here('data', 'GWAS_results', 'scoresTwo.csv')
write.GWASpoly(data = dataTwo, trait = 'Starch_ratio', filename = outFileScoresTwo, what = 'scores', delim = ',')
outFileEffectsTwo <- here('data', 'GWAS_results', 'effectsTwo.csv')
write.GWASpoly(data = dataTwo, trait = 'Starch_ratio', filename = outFileEffectsTwo, what = 'effects', delim = ',')
# Diplo-general model
outFileScoresDipGen <- here('data', 'GWAS_results', 'scoresDipGen.csv')
write.GWASpoly(data = dataDipGen, trait = 'Starch_ratio', filename = outFileScoresDipGen, what = 'scores', delim = ',')
# Diplo-general model
outFileScoresDipAdd <- here('data', 'GWAS_results', 'scoresDipAdd.csv')
write.GWASpoly(data = dataDipAdd, trait = 'Starch_ratio', filename = outFileScoresDipAdd, what = 'scores', delim = ',')
outFileEffectsDipAdd <- here('data', 'GWAS_results', 'effectsDipAdd.csv')
write.GWASpoly(data = dataDipAdd, trait = 'Starch_ratio', filename = outFileEffectsDipAdd, what = 'effects', delim = ',')
# General model
outFileScoresGen <- here('data', 'GWAS_results', 'scoresGen.csv')
write.GWASpoly(data = dataGen, trait = 'Starch_ratio', filename = outFileScoresGen, what = 'scores', delim = ',')