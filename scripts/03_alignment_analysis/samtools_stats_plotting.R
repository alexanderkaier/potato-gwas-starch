###############################
#                             #
# Importing necessary modules #
#                             #
###############################

library(tidyverse) # the god module :)
library(here)


######################
#                    #
# Importing the data #
#                    #
######################

# Loading stats data from the pure analysis without considering duplicates
statFrame <- read_csv(here("analysis", "alignment_data_GWAS", "summary_table.csv"))
# Loading stats data from the marked duplicates
statFrameMarkDup <- read_csv(here("analysis", "alignment_data_markeddup", "summary_table_markeddup.csv"))
# Loading stats data from the removed duplicates
statFrameRemDup <- read_csv(here("analysis", "alignment_data_markeddup", "summary_table_removeddup.csv"))

# Create a boxplot with 
#ggplot(data = statFrame, aes(x=factor(aligner), y=percentage_of_properly_paired_reads)) + 
#  geom_boxplot()


##########################################
#                                        #
# Extracting some basic stats for the MA #
#                                        #
##########################################

# Genome length (chr01-12)
genLength12 <- 725017384
genLengthAllChr <- 884108296

# Total number of input sequences
sum(statFrame$sequences); mean(statFrame$sequences)
totSeqs <- sum(statFrame$sequences)
# Total number and ratio of mapped reads
sum(statFrame$reads_mapped); sum(statFrame$reads_mapped)/totSeqs
# Total number and ratio of paired reads
sum(statFrame$reads_paired); sum(statFrame$reads_paired)/totSeqs
# Total number and ratio of mapped and paired reads
sum(statFrame$reads_mapped_and_paired); sum(statFrame$reads_mapped_and_paired)/totSeqs; sum(statFrame$reads_mapped) - sum(statFrame$reads_mapped_and_paired)
# Total number and ratio of properly-paired reads
sum(statFrame$reads_properly_paired); sum(statFrame$reads_properly_paired)/totSeqs
# Average quality