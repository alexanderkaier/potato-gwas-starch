#!/usr/bin/python3
# -*- coding: utf-8 -*-

####!/home/akaier/Environments/Bioinformatics/bin/python3



#####################
#                   #
# Importing modules #
#                   #
#####################

import argparse
import os
import sys
import vcf
import operator
import math
import pathlib
from collections import defaultdict, namedtuple


################################################################
#                                                              #
# Defining a new parser class and all arguments of the program #
#                                                              #
################################################################

# Define a new parser class for handling from flags or no input argument at all
class MyArgumentsParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

# Create the parser
my_parser = MyArgumentsParser(prog='vcfFilter',
                                    usage='%(prog)s [options] <vcf file>',
                                    description='Filter a given variant call file in VCF format for certain parameters',
                                    allow_abbrev=False,
                                    epilog='What do baristas wear nowadays?.. A coughy filter!')
# Specifying the version number of the filter program
my_parser.version = '1.0'

# Add the arguments
# Add the positional arguments (only one here, the input vcf file)
my_parser.add_argument('vcf',
                       metavar='<vcf file>',
                       type=str,
                       nargs='?',
                       help='the input vcf file (obligate)')

# Add the optional arguments
my_parser.add_argument('-v',
                       '--version',
                       action='version')
my_parser.add_argument('-q',
                       '--qual',
                       metavar='',
                       action='store',
                       type=float,
                       help='filtering records according to their QUAL value. The input must be an integer or a float specified \
                             to a precision up to one decimal position. Higher precision gets truncated to one decimal position')
my_parser.add_argument('-i',
                       '--info',
                       metavar='',
                       action='store',
                       type=str,
                       help='filtering records by values in the INFO area of the input file')
my_parser.add_argument('-a',
                       '--all-alleles',
                       action='store_true',
                       help='Determining, whether all alternative alleles must satisfy the criteria \
                             specified in the --info argument (concerns only multiallelic variations). \
                             If specified, records containing alleles that don´t suffice the criteria get \
                             filtered out, even if some alleles might suffice the criteria. If not specified, \
                             a record is kept, if at least one allele fulfills the criteria')
my_parser.add_argument('-f',
                       '--format',
                       metavar='',
                       action='store',
                       type=str,
                       help='filtering records by values in the FORMAT area of the input file')
my_parser.add_argument('-s',
                       '--whole-sample',
                       action='store_true',
                       help="Determining, whether all alternative alleles within a sample genotype must satisfy the criteria \
                             specified in the --format argument (concerns only multiallelic variations). \
                             If specified, samples containing alleles that don´t suffice the criteria get \
                             filtered out and replaced by a '.' in the output file, even if some alleles might suffice the criteria. If not specified, \
                             a sample is kept, if at least one allele fulfills the criteria")
my_parser.add_argument('-m',
                       '--missing-data',
                       metavar='',
                       action='store',
                       type=float,
                       help='percentage threshold of missing data for a record to be filtered \
                            (Missing data refers to number of samples not possessing a called genotype and all records \
                            possessing uncalled genotypes below the threshold are filtered out). \
                            (Specify percentage up to one decimal place (e.g. 12.5). Specified values with higher \
                            precision get truncated to one decimal position')
my_parser.add_argument('-b',
                       '--format-before-missing',
                       action='store_true',
                       help='Determing, whether missing percentage of a record or samples should be filtered first. If specified, \
                             samples not satisfying the filter criteria get filtered out BEFORE the percentage of uncalled samples of a variant is calculated. \
                             If not specified, it is the other way around, first filtering the variants based on the percentage of uncalled samples \
                             and subsequently the filtering of samples within the remaining variant. The latter approach is more native, but certain tasks \
                             might demand the former approach')
my_parser.add_argument('-o',
                       '--out',
                       metavar='',
                       action='store',
                       type=str,
                       nargs='?',
                       help='The name of the output file (path relative to execution directory of this program). \
                             This is mandatory')



#######################################
#                                     #
# Defining further necessary features #
#                                     #
#######################################

# Execute the parse_arg() method
args = my_parser.parse_args()

# Defining the necessary operators for parsing boolean operators passed as strings
ops = {
    '+' : operator.add,
    '-' : operator.sub,
    '*' : operator.mul,
    '/' : operator.truediv,
    '%' : operator.mod,
    '^' : operator.xor,
    '>': operator.gt,
    '<': operator.lt,
    '>=': operator.ge,
    '<=': operator.le,
    '=': operator.eq,
    '!=': operator.ne
}

# Defining a truncation function to truncate floats at any desired decimal place (e.g. specified float values with ridiculously high precision)
def truncate(number, decimals=0):
    '''
    Returns a value truncated to a specific number of decimal places.
    '''
    if not isinstance(decimals, int):
        raise TypeError('decimal places must be an integer.')
    elif decimals < 0:
        raise ValueError('decimal places has to be 0 or more.')
    elif decimals == 0:
        return math.trunc(number)

    factor = 10.0 ** decimals
    return math.trunc(number * factor) / factor



#############################################################################################################
#                                                                                                           #
# Error detection section to catch erroneous input by the user, print an error message and exit the program #
#                                                                                                           #
#############################################################################################################


# Checking, which flags were specified, performing error checking for those and returning useful data structures for filtering
def vcf_path_checker(relPathString):
    vcfPathString = pathlib.Path(relPathString)
    try:
        abs_path = vcfPathString.resolve(strict=True)
        return abs_path
    except FileNotFoundError as FE:
        raise FileNotFoundError("The specified input file path does not exist: '{}'".format(os.path.abspath(vcfPathString))) from FE


def qual_checker(qual_threshold):
    if qual_threshold < 0.0:
        print("Specified value for missing data is negative but must be positive. Specified value: {}%".format(qual_threshold))
        sys.exit(1)
    else:
        qual_threshold = truncate(qual_threshold, decimals=1)
    return qual_threshold


# Checking the validity of input option values for the FORMAT field
# The FORMAT field can be accessed via the vcfFile.formats attribute.
# It is stored as an ordered dictionary

# Checking for erroneous FORMAT fields in the input string that don't exist in the VCF file and for the correct number and type of passed values for each field
def format_error_checker(format_string):
    '''
    This function checks for wrongly specified fields or values for filtering samples within records concerning their
    FORMAT fields. The errors are stored in a list or dictionary which are used for printing an error message
    before exiting the program. The following misuse is catched by the funtion:
        - Specified, but not existent field(s) in the input VCF file
        - Wrong number of values for a field to be filtered
        - Wrong type of specified value for a field to be filtered
    Example: Concering a record with 1 alternative allele
    vcfFilter.py -i "GQ > 20, AO >= 10" --> no error reported, since both 'AB' and 'SAR' are valid info fields and contain each a single value in a bi-allelic variant
    vcfFilter.py -i "SART > 13, AB <= 0.4" --> error reported for specifying a non-existent field within the VCF file
    vcfFilter.py -i "SAR > 13 14, AB <= 0.4" --> error reported for specifying to many inyput values for the 'SAR' field
    '''
    format_list = format_string.split(',')
    format_dict = defaultdict()
    for condition in format_list:
        condition = condition.strip() # Strip leading and tailing whitespaces from each condition
        condition_list = condition.split(' ') # Split each condition at its spaces into three variables: 'FIELD', 'OPERATOR', 'THRESHOLD_VALUE'
        format_dict[condition_list[0]] = condition_list[1:] # The key is the format field, the two values represent the operator and the threshold provided by the user

    # Initializing the used data structures to record erroneous input
    error_list_keys = []
    error_list_fracs = []
    error_dict_values = defaultdict()
    error_dict_value_types = defaultdict()
    error_denom_dict = defaultdict()
    # Looping over all specified fields to check for erroneous input
    for key, value in format_dict.items():
        # If the input FORMAT field does not contain a division of two FORMAT fields, validate each keys input specs
        if '/' not in key:
            # If the input FORMAT field is not in the VCF file, it is added to the error_list_keys list
            if key not in vcfFile.formats.keys():
                error_list_keys.append(key)
            # If the number of values for a key is not 1 for that field in the VCF file, it is added to the error_dict_values dictionary
            if len(value) != 2: # The length of the value must be reduced by 1 to correct for the operator within the list
                error_dict_values[key] = [len(value)-1] # Saving the wrongly saved number of values for that key
            # Checking for correct input value type (distinguish between integer, float, and everything else (i.e. non-numeric))
            # If a value has not the expected type, its type is added to the error_dict_value_types dictionary, together with the expected type, for error message purposes
            if vcfFile.formats[key].type == 'Integer':
                try:
                    value[1] = int(value[1])
                except ValueError:
                    try:
                        x = float(value[1])
                        error_dict_value_types[key] = [x, vcfFile.formats[key].type]
                    except:
                        error_dict_value_types[key] = [value[1], vcfFile.formats[key].type]
            elif vcfFile.formats[key].type == 'Float':
                try:
                    value[1] = float(value[1])
                except ValueError:
                    try:
                        x = int(value[1])
                        error_dict_value_types[key] = [x, vcfFile.formats[key].type]
                    except:
                        error_dict_value_types[key] = [value[1], vcfFile.formats[key].type]

        # If the input FORMAT field is a division, validate both, numerator and denominator, and the specified value
        else:
            keystring = key.split('/')
            error_list_fracs = []
            # Checking, whether a part of the fraction is not part of the FORMAT fields in the input file and appending the error dictionary, if necessary
            for part in keystring:
                if part not in vcfFile.formats.keys():
                    error_list_fracs.append(part)
            if len(error_list_fracs) > 0:
                continue
            # The denominator of the filter field must be a single value (encoded in the VCF object as vcfFile.formats[key][2])
            elif vcfFile.formats[keystring[1]][1] != 1:
                error_denom_dict[key] = True
            # If the number of values for a key is not 1 for that field in the VCF file, it is added to the error_dict_values dictionary
            elif len(value) != 2: # The length of the value must be reduced by 1 to correct for the operator within the value list
                error_dict_values[key] = [len(value)-1] # Saving the wrongly saved number of values for that key
            # Checking for correct input value type (distinguish between integer, float, and everything else (i.e. non-numeric))
            # If a value has not the expected type, its type is added to the error_dict_value_types dictionary, together with the expected type, for error message purposes
            else:
                try:
                    value[1] = float(value[1])
                    value[1] = truncate(value[1], decimals=1)
                except ValueError:
                    try:
                        x = int(value[1])
                        error_dict_value_types[key] = [x]
                    except:
                        error_dict_value_types[key] = [value[1]]

    # Printing error messages in case of wrong input and exiting the program
    # If no error is detected, return the created info dictionary
    if len(error_list_keys) > 0:
        vcf_format_key_list = list(vcfFile.formats.keys()) # Getting all FORMAT keys of the input file for error printing
        vcf_format_key_string = ', '.join(vcf_format_key_list)
        error_key_string = ', '.join(error_list_keys) # Getting all erroneous keys of the input field for error printing
        print("Following field(s) do(es) not exist in in the FORMAT area of the input file: {}. Please check your input for the -f (--format) argument".format(error_key_string), \
              "Please refer to one of the following, available keywords or ratios thereof using a '/': {}".format(vcf_format_key_string))
        sys.exit(1)
    if len(error_list_fracs) > 0:
        vcf_format_key_list = list(vcfFile.formats.keys()) # Getting all FORMAT keys of the input file for error printing
        vcf_format_key_string = ', '.join(vcf_format_key_list)
        error_key_string = ', '.join(error_list_fracs) # Getting all erroneous keys of the input field for error printing
        print("Following field(s) do(es) not exist in in the FORMAT area of the input file: {}. Please check your input for the -f (--format) argument".format(error_key_string), \
              "Please refer to one of the following, available keywords or ratios thereof using a '/': {}".format(vcf_format_key_string))
        sys.exit(1)
    if len(error_dict_values) > 0:
        for key, values in error_dict_values.items():
            print("Number of specified values for field '{}' is {}. Expected number is {}. Please check your input for the -i (--format) argument".format(key, values[0], "1"))
        sys.exit(1)
    if len(error_dict_value_types) > 0:
        # Creating type dictionary for unambigious error messages
        type_dict = {
            int : 'Integer',
            float : 'Float',
            str : 'String'
        }
        for key, value in error_dict_value_types.items():
            if len(value) == 2:
                print("Specified value and its type for field '{}' are '{}' and '{}', respectively. Expected type is {}. Please check your input for the -i (--format) argument ".format(key, value[0], type_dict[type(value[0])], value[1]), \
                      "and the VCF documentation for further information")
            else:
                print("Specified value and its type for field '{}' are '{}' and '{}', respectively. Expected type is Float. Please check your input for the -i (--format) argument ".format(key, value[0], type_dict[type(value[0])]), \
                      "and the VCF documentation for further information")
        sys.exit(1)
    if len(error_denom_dict) > 0:
        error_denom_list = list(error_denom_dict.keys()) # Getting all FORMAT keys of the input file for error printing
        error_denom_string = ', '.join(error_denom_list)
        # Creating a list of possible denominators for error messages and troubleshoting
        denom_file_list = []
        for key, values in vcfFile.formats.items():
            if values[1] == 1:
                denom_file_list.append(key)
        print("The following field(s) were erroneously specified: {}. The denominator of a fraction can only be of length one.".format(error_denom_string), \
              "Please refer to one of the following keys for consideration as denominator in a fraction when specifying filter conditions: {}".format(denom_file_list))
        sys.exit(1)
    return format_dict


# Checking the validity of input option values for the INFO field
# The INFO field can be accessed via the vcfFile.infos attribute.
# It is stored as an ordered dictionary

# Checking for erroneous INFO fields in the input string that don't exist in the VCF file and for the correct number and type of passed values for each field
def info_error_checker(info_string):
    '''
    This function checks for wrongly specified fields or values for filtering records concerning their
    INFO fields. The errors are stored in a list or dictionary which are used for printing an error message
    before exiting the program. The following misuse is catched by the funtion:
        - Specified, but not existent field(s) in the input VCF file
        - Wrong number of values for a field to be filtered
        - Wrong type of specified value for a field to be filtered
    Example: Concering a record with 1 alternative allele
    vcfFilter.py -i "SAR > 13, AB <= 0.4" --> no error reported, since both 'AB' and 'SAR' are valid info fields and contain each a single value in a bi-allelic variant
    vcfFilter.py -i "SART > 13, AB <= 0.4" --> error reported for specifying a non-existent field within the VCF file
    vcfFilter.py -i "SAR > 13 14, AB <= 0.4" --> error reported for specifying to many inyput values for the 'SAR' field
    '''

    # Tranlate the info string into a list
    info_list = info_string.split(',')
    # Creating a dictionary from the input list
    info_dict = defaultdict()
    for condition in info_list:
        condition = condition.strip() # Strip leading and tailing whitespaces from each condition
        condition_list = condition.split(' ') # Split each condition at its spaces into three variables: 'FIELD', 'OPERATOR', 'THRESHOLD_VALUE'
        info_dict[condition_list[0]] = condition_list[1:] # The key is the info field, the two values represent the operator and the threshold provided by the user

    # Initializing the used data structures to record erroneous input
    error_list_keys = []
    error_dict_values = defaultdict()
    error_dict_value_types = defaultdict()
    # Looping over all specified fields to check for erroneous input
    for key, value in info_dict.items():
        # If the input INFO field is not in the VCF file, it is added to the error_list_keys list
        if key not in vcfFile.infos.keys():
            error_list_keys.append(key)
        # If the number of values for a key is not 1 for that field in the VCF file, it is added to the error_dict_values dictionary
        elif len(value) != 2: # The length of the value must be reduced by 2 to correct for the operator within the list
            error_dict_values[key] = [len(values)-1] # Saving both, the specified number of values and the corresponding key for error messages
        # Checking for correct input value type (distinguish between integer, float, and everything else (i.e. non-numeric))
        # If a value has not the expected type, its type is added to the error_dict_value_types dictionary, together with the expected type, for error message purposes
        elif vcfFile.infos[key][2] =='Integer':
            try:
                value[1] = int(value[1])
            except ValueError:
                error_dict_value_types[key] = [value[1], vcfFile.infos[key][2]]
        elif vcfFile.infos[key][2] =='Float':
            try:
                value[1] = float(value[1])
            except ValueError:
                error_dict_value_types[key] = [value[1], vcfFile.infos[key][2]]

    # Printing error messages in case of wrong input and exiting the program
    # If no error is detected, return the created info dictionary
    if len(error_list_keys) > 0:
        key_string = ', '.join(error_list_keys)
        print("Following field(s) do(es) not exist in in the INFO area of the input file: {}. Please check your input for the -i (--info) argument".format(key_string))
        sys.exit(1)
    elif len(error_dict_values) > 0:
        for key, values in error_dict_values.items():
            print("Number of specified values for field '{}' is {}. Expected number is {}. Please check your input for the -i (--info) argument".format(key, values[0], "1"))
        sys.exit(1)
    elif len(error_dict_value_types) > 0:
        for key, value in error_dict_value_types.items():
            print("Specified value for field '{}' is '{}'. Expected type is {}. Please check your input for the -i (--info) argument ".format(key, value[0], value[1]), \
                  "and the VCF documentation for further information")
        sys.exit(1)
    elif 'TYPE' in info_dict.keys():
        if info_dict['TYPE'][1] not in ['snp', 'mnp', 'ins', 'del', 'complex']:
            print("The specified criterion for the 'TYPE' field is '{}' and therefore not one of the following: ".format(info_dict['TYPE'][1]),\
                  "'snp', 'mnp', 'ins', 'del', 'complex'. Please check your input for the -i (--info) argument")
            sys.exit(1)
        elif info_dict['TYPE'][0] not in ['=', '!=']:
            print("The specified operator for the 'TYPE' field is '{}' and therefore not one of the following: ".format(info_dict['TYPE'][0]),\
                  "'=' or '!=', i.e. operators checking for equality. Please check your input for the -i (--info) argument")
            sys.exit(1)
    else:
        return info_dict


# Parsing and checking for the validity of the -m input
def missing_data_checker(missing_const):
    if missing_const < 0.0 or missing_const > 100.0:
        print("Specified value for missing data is not within the expected range of 0-100%. Specified value: {}%".format(missing_const))
        sys.exit(1)
    else:
        missing_const = truncate(missing_const, decimals=1)
    return missing_const


# Check, if specified output path exists and creating new path if not
def out_path_checker(path_string):
    #vcfOutString = os.path.abspath(args.out)
    vcfOutAbsString = os.path.abspath(args.out)
    # Create a list for seperating specified path into directory and filename
    vcfOutList = vcfOutAbsString.split('/')
    vcfOutDir = '/'.join(vcfOutList[:-1])
    vcfOutFileName = vcfOutList[-1]
    # Check if the specified output directory already exists and create it if not
    if not os.path.exists(vcfOutDir):
        os.makedirs(vcfOutDir)
        open(vcfOutAbsString, 'w').close()
    # If the output file doesn't exist yet, it will be created
    elif not os.path.isfile(vcfOutAbsString):
        open(vcfOutAbsString, 'w').close()
    # If the output file already exists, the user has to decide whether to overwrite it or to terminate the program
    else:
        valAns = ['y', 'n']
        ans = input("Output file already exists. Do you want to overwrite it [y/n]? ")
        while ans not in valAns:
            ans = input("Your answer is invalid. Please type 'y', or 'n' for overwriting the output file and continuing or keeping the output file and terminating the program, respectively! ")
        if ans == 'n':
            print("Output file kept and program terminated. Have a good day!")
            sys.exit(0)
        else:
            open(vcfOutAbsString, 'w').close() # Overwriting existent file if desired
    return vcfOutAbsString


# Check, if no input argument was given and print help message, if so
if len(sys.argv) == 1:
    my_parser.print_help(sys.stderr)
    sys.exit(1)



###################################################################################################################################################
#                                                                                                                                                 #
# Checking for specified arguments, performing error checking via pre-defined functions and saving specified arguments in variables for filtering #
#                                                                                                                                                 #
###################################################################################################################################################

# Saving the input file string to a variable, if specified
if args.vcf is not None:
    vcfPath = vcf_path_checker(args.vcf)
    vcfString = str(vcfPath)
    vcfFile = vcf.Reader(filename=vcfString)
else:
    my_parser.print_help(sys.stderr)
    sys.exit(1)
# Saving the quality specification as a variable, if defined
if args.qual is not None:
    qual_threshold = qual_checker(args.qual)
# Saving the format specifications as a variable, if defined
if args.format is not None:
    format_dict = format_error_checker(args.format)
# Saving the info specifications as a variable, if defined
if args.info is not None:
    info_dict = info_error_checker(args.info)
# Saving the specifications concerning tolerable missing data as a variable, if defined
if args.missing_data is not None:
    missing_perc = missing_data_checker(args.missing_data)
# Saving the output path as a variable, if defined
if args.out is not None:
    vcfOutPath = out_path_checker(args.out)



###########################################################################################################
#                                                                                                         #
# Defining function for filtering VCF file variants for later application in a specified, efficient order #
#                                                                                                         #
###########################################################################################################

# Setting the first record in the file preliminary for checking the FORMAT specifications
rec = next(vcfFile)
# Creating key list of the the FORMAT area to check for actual FORMAT fields of the records
format_key_list = rec.FORMAT.split(':')
# Initializing error list for record FORMAT areas, that are specified, but not present in the input instance file
rec_format_error_list = []
# Appending the error list by every key, that was specified and exists theoretically in the VCF FORMAT, but was not implemented in the input file (e.g. 'GL' or 'MIN_DP' in this case)
if 'format_dict' in vars():
    for key in format_dict.keys():
        if '/' not in key:
            if key not in format_key_list:
                rec_format_error_list.append(key)
        else:
            for part in key.split('/'):
                if part not in format_key_list:
                    rec_format_error_list.append(part)
    if len(rec_format_error_list) > 0:
        rec_format_error_string = ', '.join(rec_format_error_list)
        print("Following specified field(s) were not generated by the variant calling program, but can be created by specifying the according argument during variant call: {}".format(rec_format_error_string))
        sys.exit(1)



###########################################################################################################
#                                                                                                         #
# Defining function for filtering VCF file variants for later application in a specified, efficient order #
#                                                                                                         #
###########################################################################################################


# Function for filtering records within the VCF input file according to their quality
def qual_filter(rec, qual):
    if rec.QUAL >= qual:
        passed = True
    else:
        passed = False
    return passed

# Function for filtering records within the VCF input file according to the percentage of missing data
def missing_filter(rec, missing_perc):
    if (1 - rec.call_rate) <= missing_perc:
        passed = True
    else:
        passed = False
    return passed


# Function for filtering records within the VCF input file according to the INFO fields
def info_filter(rec, info_dict):
    # Initializing a dictionary to store all keys and lists containing boolean values, whether the vcf entry or allele passed the threshold
    info_dict_checked = defaultdict()
    for key, values in info_dict.items():
        items = []
        # if this record has only one value for that field, this value will be compared directly with the specified threshold
        if not hasattr(rec.INFO[key], '__iter__'):
            #items.append(ops[values[0]](rec.INFO[key],values[1]))
            items += len(rec.ALT) * [ops[values[0]](rec.INFO[key],values[1])]
            info_dict_checked[key] = items
        # If the record has a list of multiple values at a position(one for each allele) then all of them get compared with the specified threshold
        # and the result saved as a list value for that field in a dictionary
        else:
            for allele in range(len(rec.INFO[key])):
                items.append(ops[values[0]](rec.INFO[key][allele],values[1]))
                info_dict_checked[key] = items

    big_list = []
    # If the -a flag is specified, all values in the dictionary must equal True
    if args.all_alleles:
        for key, allele_list in info_dict_checked.items():
            X = all(info_dict_checked[key])
            if not X:
                passed = False
                break
            else:
                passed = True
    # If the -a flag is not set, then every allele has to be checked for satisfying ALL specified criteria
    else:
        for key, allele_list in info_dict_checked.items():
            big_list.append(allele_list)
        # Iterating over all alleles within the 'big_list' to check, if at least one passes all filter criteria
        for i in range(len(big_list[0])):
            allele_check = []
            for allele in big_list:
                allele_check.append(allele[i])
            X = all(allele_check)
            if X:
                passed = True
                break
            else:
                passed = False
    return passed


# Function for filtering records within the VCF input file according to the FORMAT fields
def format_filter(rec, format_dict):
    # Creating a list of all keys in the FORMAT area of the record
    key_list = rec.FORMAT.split(':')
    # Looping over each sample in the record
    for sample in rec.samples:
        # Pass the sample, if its genotype was not called at that position
        if not sample.called:
            continue
        else:
            # Creating a boolean dictionary with the keys being the filtered FORMAT fields and the values attaining True or False for each alternative allele and key
            format_dict_checked = defaultdict()
            for key, values in format_dict.items():
                # Initializing an empty list for the field (its length will be the number of alternative alleles of that record)
                items = []
                if '/' not in key:
                    # if this record has only one value for that field, this value will be compared directly with the specified threshold
                    if not hasattr(sample[key], '__iter__'):
                        # If there is only one value for that field (i.e. 'DP'), its value is multiplied by the number of alternative alleles and saved as value in the dictionary (this simplifies the subsequent filtering)
                        items += len(rec.ALT) * [ops[values[0]](sample[key],values[1])]
                        format_dict_checked[key] = items
                    else:
                        for allele in range(len(sample[key])):
                            items.append(ops[values[0]](sample[key][allele],values[1]))
                            format_dict_checked[key] = items
                # Comparing specified ratios with calculated ones
                else:
                    keystring = key.split('/')
                    # Check if the numerator consists only of a single value in the input file
                    if not hasattr(sample[keystring[0]], '__iter__'):
                        # If there is only one value for the numerators field (i.e. 'DP'), the desired ratio is calculated and compared with the specified threshold. The boolean value gets multiplied by the number of alternative alleles called at that record
                        numer = float(sample[keystring[0]])
                        denom = float(sample[keystring[1]])
                        ratio_perc = (numer / denom) * 100
                        items += len(rec.ALT) * [ops[values[0]](ratio_perc, values[1])]
                        format_dict_checked[key] = items
                    else:
                        for allele in range(len(rec.ALT)):
                            numer = float(sample[keystring[0]][allele])
                            denom = float(sample[keystring[1]])
                            ratio_perc = (numer / denom) * 100
                            items.append(ops[values[0]](ratio_perc, values[1]))
                            format_dict_checked[key] = items


        # If the -s flag is specified, all values in the dictionary for a single sample must equal True
        if args.whole_sample:
            for key, allele_list in format_dict_checked.items():
                X = all(allele_list)
                # If not all values are true, the genotype ('GT') of the sample is changed to '.' and every other field is changed to None
                # Otherwise, nothing happens to the samples genotype information
                if not X:
                    # Turn the named tuple of the sample namespace (sample.data) into a dictionary for changing values and replace the old named tuple with it
                    call_tup = sample.data
                    call_dict = call_tup._asdict()
                    # Looping over all FORMAT keys in the record and updating the CallData dictionary
                    for field in rec.FORMAT.split(':'):
                        if field == 'GT':
                            call_dict['GT'] = '.'
                        else:
                            call_dict[field] = None
                    # Transforming the 'blank' dictionary back into a named tuple called 'CallData' and assigning it to sample.data, thereby overwriting the previous values for sample.data
                    data_tup = namedtuple('CallData', call_dict)
                    new_call_tup = data_tup(**call_dict)
                    # Assigning new values to genotype and whether it was called for potential subsequent filtering
                    sample.data = new_call_tup
                    sample.called = False
                    break
                else:
                    continue
        # If the -s flag is not set, then every allele of a sample has to be checked for satisfying ALL specified criteria, but one is enough to pass the filter
        else:
            big_list = []
            for key, allele_list in format_dict_checked.items():
                big_list.append(allele_list)
            # Iterating over all alleles within the 'big_list' to check, if at least one passes all filter criteria
            for i in range(len(big_list[0])):
                allele_check = []
                for allele in big_list:
                    allele_check.append(allele[i])
                X = all(allele_check)
                if X:
                    passed = True
                    break
                else:
                    passed = False
            # If no allele passed the criteria, the sample gets rendered uncalled
            if not passed:
                # Turn the named tuple of the sample namespace (sample.data) into a dictionary for changing values and replace the old named tuple with it
                call_tup = sample.data
                call_dict = call_tup._asdict()
                # Looping over all FORMAT keys in the record and updating the CallData dictionary
                for field in rec.FORMAT.split(':'):
                    if field == 'GT':
                        call_dict['GT'] = '.'
                    else:
                        call_dict[field] = None
                # Transforming the 'blank' dictionary back into a named tuple called 'CallData' and assigning it to sample.data, thereby overwriting the previous values for sample.data
                data_tup = namedtuple('CallData', call_dict)
                new_call_tup = data_tup(**call_dict)
                # Assigning new values to genotype and whether it was called for potential subsequent filtering
                sample.data = new_call_tup
                sample.called = False




#######################################################################################
#                                                                                     #
# Actual filering based on validation of input and filter functions implemented above #
#                                                                                     #
#######################################################################################



# Creating the output object using the output string
vcfOutFile = vcf.Writer(open(vcfOutPath, 'w'), vcfFile)

# Looping over all records in the VCF input file, thereby passing one filter after another and directly skipping the current variant, if the variable 'passed' attains the value 'False' for any filter
for rec in vcfFile:
    if 'qual_threshold' in vars():
        passed = qual_filter(rec=rec, qual=qual_threshold)
        if not passed:
            continue
    elif 'info_dict' in vars():
        passed = info_filter(rec=rec, info_dict=info_dict)
        if not passed:
            continue
    elif args.format_before_missing:
        if 'format_dict' in vars():
            format_filter(rec=rec, format_dict=format_dict)
        elif 'missing_perc' in vars():
            passed = missing_filter(rec=rec, missing_perc=missing_perc)
            if not passed:
                continue
    else:
        if 'missing_perc' in vars():
            passed = missing_filter(rec=rec, missing_perc=missing_perc)
            if not passed:
                continue
        elif 'format_dict' in vars():
            format_filter(rec=rec, format_dict=format_dict)

    # Writing selected variant either to stdout or output file, if specified
    vcfOutFile.write_record(rec)
print('Done!')

