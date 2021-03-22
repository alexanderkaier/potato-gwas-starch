#!/bin/bash

# This code copies the first 50 sample folders from the 192 samples from the GWAS folder into this project
# in order to establish the benchmark using those 50 samples exemplarily.
find ../../GWAS/data/GWAS_data/Final_BIG_GBS/RE_processed/* -maxdepth 0 | sort -k1 | head -n50 | xargs cp -r -t ../data/initial_seqs