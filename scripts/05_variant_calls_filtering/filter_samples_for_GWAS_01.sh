

# Filtering SNP VCF file for samples that were phenotyped
# bcftools view -S sample

# In- and output for the program
INPUT=../../data/Phenotype_data/GWAS_input.csv
INTERM_SAMPLE_FILE=samples.txt
VCF_IN_FILE=../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs.vcf
VCF_OUT_FILE=../../data/variant_calls_GWAS_filtered/freebayes_192_samples_chr01-chr12_SNPs_135-samples.vcf
IN_PATH_SUBSET=../../data/variant_calls_GWAS_filtered

# Creating the intermediate samples TXT file from all samples in the phenotype GWAS input file
{
    read REPLY
    while IFS=, read -r individual trait
    do
        echo "$individual" >> $INTERM_SAMPLE_FILE
    done
} < $INPUT


# Actual filtering of samples
bcftools view --samples-file samples.txt $VCF_IN_FILE > $VCF_OUT_FILE

# Deleting the intermediate samples TXT file
rm $INTERM_SAMPLE_FILE

# Subsetting approx. 100,000 SNPs randomly from the SNPs file with 135 samples
bcftools view $VCF_OUT_FILE | vcfrandomsample --rate 0.054 \
    > $IN_PATH_SUBSET/freebayes_192_samples_chr01-chr12_SNPs_135-samples_subset100000.vcf