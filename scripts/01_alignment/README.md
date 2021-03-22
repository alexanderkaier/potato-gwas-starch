# *Alignment of short reads to a reference genome*

This folder contains all scripts necessary for indexing reference fasta files as well as alignment of trimmed reads using bwa mem and Bowtie2.<br/>No directories need to be initialized, as every script does it on its own when called.

The basic workflow is usually the same, regardless of the used tool. First, the reference FASTA file needs to be indexed only ONCE.<br/>Every subsequent alignment can be done using the reference index generated initially as reference for alignment.

## Workflow using bwa mem

From the directory containing this README.md file, run the commands in the following order:
- ./bwa-mem-reference-indexing.sh
- ./bwa-mem-alignment.sh

