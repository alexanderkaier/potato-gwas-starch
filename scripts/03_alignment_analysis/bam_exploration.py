#!/usr/bin/python

# Importing used libraries.
import pysam
import os
from pprint import pprint

# Defining the input BAM file directory
twd = os.path.abspath('../../data/alignment_data/bwa-mem/Sample_P1-A01-001-AD/Sample_P1-A01-001-AD_sorted.bam')
print(twd)

# creating bam object and extracting headers
#bam = pysam.AlignmentFile(twd)
bam = pysam.AlignmentFile(twd)
headers = bam.header

# printing header attributes
for record_type, records in headers.items():
	print(record_type)
	for i, record in enumerate(records):
		print('\t%d' % (i+1))
		if type(record) == str:
			print('\t\t%s' % record)
		elif type(record) == dict:
			for field, value in record.items():
				print('\t\t%s' % (field, value))

# Cecking a single record
for rec in bam:
	if rec.cigarstring.find('M') > -1 and \
	   rec.cigarstring.find('S') > -1 and \
	   not.rec.is_unmapped and \
	   not rec.mate_is_unmapped:
	   break
print(rec.query_name, rec.referece_id, bam.getrname(rec.referece_id), rec.reference_start, rec.reference_end)
print(rec.cigarstring)
print(rec.query_alignment_start, rec.query_alignment_end, rec.query_alignment_length)
print(rec.next_reference_id, rec.next_reference_start, rec.template_length)
print(rec.is_paired)