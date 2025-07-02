#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/11/28
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script computes TPM from mapping and count of the reads on the 
individual CDS. 
'''
# ---------------------------------------------------------------------------
import os
import argparse
import pandas as pd
import logging
logging.basicConfig(level=logging.INFO)
# ---------------------------------------------------------------------------
# Definitions

# ---------------------------------------------------------------------------

# =============
# GET ARGUMENTS
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = '''
This script computes TPM from mapping and count of the reads on the 
individual CDS. 
'''
)

parser.add_argument("-i", "--input", help="Tsv files containing 3 column: Gene, Length, NbMappedReads. The files must be named [Strain].ReadCounts.tsv. ", required=True, nargs='+')
parser.add_argument("-o", "--output", help="Name of the output file tsv file (matrix of TPM)", required=True)

# Read arguments from the command line
args = parser.parse_args()
inputPaths = [os.path.abspath(x) for x in args.input]
outputPath = os.path.abspath(args.output)

output = []
for path in inputPaths:
	# Get strain name
	strain = path.split('/')[-1][0:-15]
	logging.info(f'Processing file: {path.split("/")[-1]}, strain: {strain}')
	# Read input file
	data = pd.read_csv(path, sep = '\t')
	data = data.loc[data['Gene'] != '*',:]
	# Compute TPM
	data['RPK'] = data['NbMappedReads'] / data['Length']
	data['SumRPK'] = sum(data['RPK'])
	data['TPM'] = data['RPK'] / (data['SumRPK'] / 1000000)
	# Format dataframe
	data.index = data['Gene']
	data = data.loc[:,['TPM']].rename(columns = {'TPM': strain})
	output += [data]

# Write output file
output = pd.concat(output, axis = 1)
output.to_csv(outputPath, sep = '\t')

