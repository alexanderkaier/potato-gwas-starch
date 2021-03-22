library(ldsep)
library(here)
library(qqplotr)


####################################
#                                  #
# Reading the highest scoring SNPs #
#                                  #
####################################

# General path to best scoring SNP and their LD linked SNPs files
AddSNPs <- as.character(unlist(read.csv(here('analysis', 'GWAS_results', 'LDSNPsAdd.csv'), header = F), use.names = F))
OneSNPs <- as.character(unlist(read.csv(here('analysis', 'GWAS_results', 'LDSNPsOne.csv'), header = F), use.names = F))
TwoSNPs <- as.character(unlist(read.csv(here('analysis', 'GWAS_results', 'LDSNPsTwo.csv'), header = F), use.names = F))
DipGenSNPs <- as.character(unlist(read.csv(here('analysis', 'GWAS_results', 'LDSNPsDipGen.csv'), header = F), use.names = F))
DipAddSNPs <- as.character(unlist(read.csv(here('analysis', 'GWAS_results', 'LDSNPsDipAdd.csv'), header = F), use.names = F))
GenSNPs <- as.character(unlist(read.csv(here('analysis', 'GWAS_results', 'LDSNPsGen.csv'), header = F), use.names = F))


#######################################
#                                     #
# Creating manhattan plots with qqman #
#                                     #
#######################################

# Create dataframes from the output of GWASpoly for each gene action model
AddFrame <- cbind(dataAdd@map$Marker, dataAdd@map$Chrom, dataAdd@map$Position, dataAdd@scores$Starch_ratio)
colnames(AddFrame) <- c("SNP", "CHR", "BP", "score")
AddFrame$CHR <- substr(AddFrame$CHR, 4,5)
AddFrame$SNP <- as.character(AddFrame$SNP)
AddFrame$CHR <- as.numeric(AddFrame$CHR)

OneFrame <- cbind(dataOne@map$Marker, dataOne@map$Chrom, dataOne@map$Position, dataOne@scores$Starch_ratio)
colnames(OneFrame) <- c("SNP", "CHR", "BP", "score1", "score2")
OneFrame$CHR <- substr(OneFrame$CHR, 4,5)
OneFrame$SNP <- as.character(OneFrame$SNP)
OneFrame$CHR <- as.numeric(OneFrame$CHR)
OneFrame$score <- rep(NA, nrow(OneFrame))
for (i in 1:nrow(OneFrame)) {
  if (!is.na(OneFrame$score1)[i] & is.na(OneFrame$score2)[i]) {
    OneFrame$score[i] = OneFrame$score1[i]
  }
  else if (is.na(OneFrame$score1)[i] & !is.na(OneFrame$score2)[i]) {
    OneFrame$score[i] = OneFrame$score2[i]
  }
  else if (!is.na(OneFrame$score1)[i] & !is.na(OneFrame$score2)[i]) {
    if (OneFrame$score1[i] >= OneFrame$score2[i]) {
      OneFrame$score[i] = OneFrame$score1[i]
    }
    else {
      OneFrame$score[i] = OneFrame$score2[i]
    }
  }
  else {
    OneFrame$score[i] = NA
  }
}

TwoFrame <- cbind(dataTwo@map$Marker, dataTwo@map$Chrom, dataTwo@map$Position, dataTwo@scores$Starch_ratio)
colnames(TwoFrame) <- c("SNP", "CHR", "BP", "score1", "score2")
TwoFrame$CHR <- substr(TwoFrame$CHR, 4,5)
TwoFrame$SNP <- as.character(TwoFrame$SNP)
TwoFrame$CHR <- as.numeric(TwoFrame$CHR)
TwoFrame$score <- rep(NA, nrow(TwoFrame))
for (i in 1:nrow(TwoFrame)) {
  if (!is.na(TwoFrame$score1)[i] & is.na(TwoFrame$score2)[i]) {
    TwoFrame$score[i] = TwoFrame$score1[i]
  }
  else if (is.na(TwoFrame$score1)[i] & !is.na(TwoFrame$score2)[i]) {
    TwoFrame$score[i] = TwoFrame$score2[i]
  }
  else if (!is.na(TwoFrame$score1)[i] & !is.na(TwoFrame$score2)[i]) {
    if (TwoFrame$score1[i] >= TwoFrame$score2[i]) {
      TwoFrame$score[i] = TwoFrame$score1[i]
    }
    else {
      TwoFrame$score[i] = TwoFrame$score2[i]
    }
  }
  else {
    TwoFrame$score[i] = NA
  }
}

DipGenFrame <- cbind(dataDipGen@map$Marker, dataDipGen@map$Chrom, dataDipGen@map$Position, dataDipGen@scores$Starch_ratio)
colnames(DipGenFrame) <- c("SNP", "CHR", "BP", "score")
DipGenFrame$CHR <- substr(DipGenFrame$CHR, 4,5)
DipGenFrame$SNP <- as.character(DipGenFrame$SNP)
DipGenFrame$CHR <- as.numeric(DipGenFrame$CHR)

DipAddFrame <- cbind(dataDipAdd@map$Marker, dataDipAdd@map$Chrom, dataDipAdd@map$Position, dataDipAdd@scores$Starch_ratio)
colnames(DipAddFrame) <- c("SNP", "CHR", "BP", "score")
DipAddFrame$CHR <- substr(DipAddFrame$CHR, 4,5)
DipAddFrame$SNP <- as.character(DipAddFrame$SNP)
DipAddFrame$CHR <- as.numeric(DipAddFrame$CHR)

GenFrame <- cbind(dataGen@map$Marker, dataGen@map$Chrom, dataGen@map$Position, dataGen@scores$Starch_ratio)
colnames(GenFrame) <- c("SNP", "CHR", "BP", "score")
GenFrame$CHR <- substr(GenFrame$CHR, 4,5)
GenFrame$SNP <- as.character(GenFrame$SNP)
GenFrame$CHR <- as.numeric(GenFrame$CHR)


png(filename=here('Figures', 'manhattan_plots_all.png'), res=130, width = 1000, height = 1300)
# Defining layout for big plot
#par(bg = 'white', mfrow = c(6,1))
layout(matrix(c(1,2,3,4,5,6), byrow = F),
       heights =c(6,6,6,6,6,7))

# Creating the manhattan plots
par(mar = c(2, 4.2, 0.5, 2))
manAdd <- manhattan(AddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08,suggestiveline = 3.25,
                    highlight = AddSNPs, ylim = c(0, 7))
#mtext(side = 2, x=-15000000,y=6, adj = c(0,1), labels = c("(a)"))

title("Additive", line = -1.5)
#text(x=300000000,y=6, adj = c(0.5, 0.8), labels = c("Additive"))
manOne <- manhattan(OneFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 4.9, suggestiveline = 3.13,
                    highlight = OneSNPs, ylim = c(0, 7))
title("Simplex", line = -1.5)
manTwo <- manhattan(TwoFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 4.79, suggestiveline = 2.91,
                    highlight = TwoSNPs, ylim = c(0, 7))
title("Duplex", line = -1.5)
manDipGen <- manhattan(DipGenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08, suggestiveline = 3.68,
                       highlight = DipGenSNPs, ylim = c(0, 7))
title("Diploid-general", line = -1.5)
manDipAdd <- manhattan(DipAddFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08, suggestiveline = 3.38,
                       highlight = DipAddSNPs, ylim = c(0, 7))
title("Diploid-additive", line = -1.5)
par(mar = c(4, 4.2, 0.5, 2))
manGen <- manhattan(GenFrame, p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08, suggestiveline = 3.9,
                    highlight = GenSNPs, ylim = c(0, 7))
title("General", line = -1.5)
dev.off()


# Creating close view of one highlighted SNP with all SNPs in LD
manhattan(subset(AddFrame, CHR == 9) , p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08, suggestiveline = 3.9,
                            highlight = AddSNPs, ylim = c(0, 6), xlim = c(3882536, 7222536))
abline(v = 4882536, lty = 2)
abline(v = 6222536, lty = 2)

manhattan(subset(OneFrame, CHR == 9) , p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 4.9, suggestiveline = 3.13,
          highlight = OneSNPs, ylim = c(0, 7), xlim = c(4e+06, 8e+06))
abline(v = 55797680, lty = 1)

##################################
#                                #
# p-values and inflation factors #
#                                #
##################################

# Calculating and plotting inflation rates for all models
pAdd <- 10^-na.omit(unlist(dataAdd@scores))
lambdaAdd <- median(qchisq(pAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pOne <- 10^-na.omit(unlist(dataOne@scores))
lambdaOne <- median(qchisq(pOne, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pTwo <- 10^-na.omit(unlist(dataTwo@scores))
lambdaTwo <- median(qchisq(pTwo, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pDipGen <- 10^-na.omit(unlist(dataDipGen@scores))
lambdaDipGen <- median(qchisq(pDipGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pDipAdd <- 10^-na.omit(unlist(dataDipAdd@scores))
lambdaDipAdd <- median(qchisq(pDipAdd, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
pGen <- 10^-na.omit(unlist(dataGen@scores))
lambdaGen <- median(qchisq(pGen, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)


####################
#                  #
# qq-plot plotting #
#                  #
####################

qqAdd <- qq.plot(dataAdd) + 
  ggtitle("Additive") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = "none",
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  annotate("text", label = "lambda[GC]^A == '1.01'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.91), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

qqOne <- qq.plot(dataOne) + 
  ggtitle("Simplex") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.5,0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        strip.text.x = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  annotate("text", label = "lambda[GC]^S == '1.00'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.9), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)

qqTwo <- qq.plot(dataTwo) + 
  ggtitle("Duplex") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
        legend.position = c(0.5,0.9),
        legend.title = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10),
        strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^D == '1.01'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=4.79), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)
 
qqDipGen <- qq.plot(dataDipGen) + 
  ggtitle("Diploid-general") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
       legend.position = "none",
       axis.title.x = element_text(size = 10),
       axis.title.y = element_text(size = 10),
       strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^DG == '0.98'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=5.08), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)
 
qqDipAdd <- qq.plot(dataDipAdd) + 
  ggtitle("Diploid-additive") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
       legend.position = "none",
       axis.title.x = element_text(size = 10),
       axis.title.y = element_text(size = 10),
       strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^DA == '1.00'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=5.08), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)
 
qqGen <- qq.plot(dataGen) + 
  ggtitle("General") + 
  theme(plot.title = element_text(hjust = 0.5, size = 15),
       legend.position = "none",
       axis.title.x = element_text(size = 10),
       axis.title.y = element_text(size = 10),
       strip.text.x = element_blank()) + 
  annotate("text", label = "lambda[GC]^G == '0.93'", parse = TRUE, x = 0.7, y = 6, size = 4.5) + 
  geom_hline(aes(yintercept=5.08), color = "#1B9E77") + 
  xlim(0, 5) + 
  ylim(0, 6.5)
 
png(filename=here('Figures', 'qq_plots_all.png'), res=130, width = 1000, height = 1300)
# Defining layout for big plot
#par(bg = 'white', mfrow = c(6,1))
#layout(matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2, byrow = T),
#       heights =c(6,6,6,6,6,7))
plot_grid(qqAdd, qqOne, qqTwo, qqDipGen, qqDipAdd, qqGen,
          nrow = 3, ncol = 2, labels = "auto", label_size = 20)
dev.off()


#########################
#                       #
# Zoomed manhattan plot #
#                       #
#########################

svg(filename=here('Figures', 'manhattan_zoom.svg'), pointsize=12)
layout(matrix(c(1,2), byrow = F),
       heights =c(6,6,6,6,6,7))
mansub9 <- manhattan(subset(GenFrame, CHR == 9) , p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08, suggestiveline = 3.9,
                    highlight = GenSNPs, ylim = c(0, 7), xlab = "General model: Chromosome 9")
abline(v = 4926242, lty = 2)
abline(v = 6266242, lty = 2)
mansubLD <- manhattan(subset(GenFrame, CHR == 9) , p = "score", logp = F, ylab = expression(paste("-log"[10],"(p)")), genomewideline = 5.08, suggestiveline = 3.9,
          highlight = GenSNPs, ylim = c(0, 7), xlim = c(3926242, 7266242), xlab = "Chromosome 9 position [bp]")
abline(v = 4926242, lty = 2)
abline(v = 6266242, lty = 2)
dev.off()
