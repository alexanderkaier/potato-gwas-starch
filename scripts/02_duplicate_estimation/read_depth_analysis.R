#########################################
#                                       #
# Analysis of depth and genome coverage #
#                                       #
#########################################

# Import libraries
library(tidyverse)
library(here)

# Loading read depth data for one sample
depthtot = read.table(here('Put path/sample_depth.txt here'), sep = '\t', header = FALSE)
colnames(depthtot) = c("Chromosome", "Position", "Depth")
depthtot$Position = as.character(depthtot$Position)

# Reading the length of the reference genome
chrLen = read.table(here("analysis/ref_genome_length.txt"))
chrLen = head(chrLen,-1)
colnames(chrLen) = c("Chromosome", "Basepairs")
refLen = sum(chrLen$V2)

# Histogram
ggplot(depthMarked %>% filter(Chromosome == "chr01", Depth >= 48), aes(x=Depth)) + 
  geom_histogram()

# Fraction of bases of the whole genome covered:
(nrow(depthtot) / refLen) * 100
(nrow(depthRemoved) / refLen) * 100

# Fraction of bases with depth >= 48
nrow(depthtot) * (nrow(depthtot %>% filter(Depth >= 48)) / nrow(depthtot))
