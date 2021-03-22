#!/bin/python


'''
Author's comment:
This script was written to reformat attribute keys of the GFF3 file that did not conform with formatting standards
'''

# Load the libraries
import pandas as pd
import pathlib
import os
import gffpandas.gffpandas as gffpd
import csv

# Creating in- and output file names
inPath = "../../data/Reference_genome_and_annotation_file/PGSC_DM_V403_genes.gff"
outPath = "../../data/Reference_genome_and_annotation_file/PGSC_DM_V403_genes_case-corrected.gff"

# Create absolute path and read the file
absInPath = os.path.abspath(inPath)
rawFile = gffpd.read_gff3(absInPath)

# Defining the function for changing the attributes column of the gff3 attributes data frame column
def attrChange(attrString):
    # Create list of official GFF3 atrtributes from https://www.ncbi.nlm.nih.gov/datasets/docs/about-ncbi-gff3/
    # This is needed to check for invalid cases in the input file
    offList = ['ID', 'Parent', 'Dbxref', 'Name', 'Note', 'Is_circular']
    attrList = attrString.split(';')
    # Check for idiotic sepearation of values within the ';'-separated list
    ''' # Commented out because 'gff3ToGenePred' can't handle multiple names per annotation
    for j in range(len(attrList)):
        if '=' not in attrList[j]:
            attrList[j-1] = attrList[j-1] + ';' + attrList[j]
    ''' 
    dummyList = attrList
    for element in dummyList:
        # Remove second names from the attribute list
        if '=' not in element:
            attrList.remove(element)
    dummyList = attrList
    for element in dummyList:
        # Remove second names from the attribute list
        if '=' not in element:
            attrList.remove(element)
        # Check for name values and replace ',' with ' '
    for i in range(len(attrList)):
        subList = attrList[i].split('=')
        #print(subList)
        if subList[0].upper() == 'NAME':
            subList[1] = subList[1].replace(',', ' ')
            # Stitching the element back together, but now without ',' in the name value field
            attrList[i] = '='.join(subList)
    for i in range(len(attrList)):
        # Extract the key from each individual annotation
        attrKey = attrList[i].split('=')[0]
        # Check, if the attributes are in the list of official ones or not. If so, uppercase the keys first letter, otherwise lowercase it
        if attrKey.upper() in (item.upper() for item in offList) and attrKey[0].islower():
            dummy = attrKey[0].upper() + attrKey[1:]
            attrList[i] = dummy + '=' + attrList[i].split('=')[1]
        elif attrKey.upper() not in (item.upper() for item in offList) and attrKey[0].isupper():
            dummy = attrKey[0].lower() + attrKey[1:]
            attrList[i] = dummy + '=' + attrList[i].split('=')[1]
    newAttrString = ';'.join(map(str, attrList))
    return newAttrString

# Apply the function to a subset of annotations for testing purposes
#rawFile.df.loc[255860:255870,:]['attributes'] = rawFile.df.loc[255860:255870,:]['attributes'].apply(lambda attrString: attrChange(attrString))
#dummyColumn = rawFile.df.loc[:10,:]['attributes'].apply(lambda attrString: attrChange(attrString))
#rawFile.df.loc[:10,:]['attributes'] = dummyColumn
rawFile.df.iloc[:,8] = rawFile.df.iloc[:,8].apply(lambda attrString: attrChange(attrString))

# Save the input file as new GFF3 file
rawFile.df.to_csv(outPath, sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)
