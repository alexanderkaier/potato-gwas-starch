# Command to print out the name of each entry (chromosome in this case) in a fasta file as well as the length of its corresponding sequence

OUT_PATH=../analysis

cat \
	../data/Reference_genome_and_annotation_file/potato_dm_v404_all_pm_un.fasta \
	| awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' \
	> ../analysis/ref_genome_length.txt
