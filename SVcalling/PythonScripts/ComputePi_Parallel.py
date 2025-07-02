#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/10/04
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script computes a Structural Pi from a SV VCF file over a sliding window. 
'''
# ---------------------------------------------------------------------------
import gzip
import pandas as pd
import argparse
import os
import logging
import itertools
from multiprocessing import Pool
logging.basicConfig(level=logging.NOTSET)
# ---------------------------------------------------------------------------
# Definitions
def harmonicNumber(n):
	k = list(range(1, n+1))
	return(sum([1/x for x in k]))

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="Path of the VCF file (gzipped or not)", required=True)
parser.add_argument("-T", "--types", help="Type of SV to keep in the VCF (any of INS DEL CONTRE DUP INV TRA, multiple allowed)", nargs = '+')
parser.add_argument("-w", "--window", help="Window size for the calculation of Pi (bp)", required=True, type = int)
parser.add_argument("-s", "--sliding", help="Length of the sliding between windows (bp)", required=True, type = int)
parser.add_argument("-o", "--output", help="Output file (VCF format, not gzipped)", required=True)
parser.add_argument("-t", "--threads", help="Number of threads on which to parallelize", type = int, default = 1)

# Read arguments from the command line
args = parser.parse_args()

vcfPath=os.path.abspath(args.vcf)
types = args.types
window = args.window
sliding = args.sliding
outputPath=os.path.abspath(args.output)
nthreads = args.threads

#vcfPath = '04_PopulationVCF/SV.1087Samples.vcf.gz'


# Read VCF =============================================================
logging.info('Read VCF file')

# Get header
header = []
if vcfPath.endswith('.gz'):
	with gzip.open(vcfPath, "rb") as vcf:
		for line in vcf:
			line = line.decode()
			if line.startswith('#CHROM'):
				colnames = line.strip().split('\t')
				break
			else:
				header += [line]
else:
	with open(vcfPath, "r") as vcf:
		for line in vcf:
			if line.startswith('#CHROM'):
				colnames = line.strip().split('\t')
				break
			else:
				header += [line]

# Change ALT colnames to ALT ALLELE to avoid confonding with ALT strain
colnames[4] = 'ALT_ALLELE'

# Read vcf
VCF = pd.read_csv(vcfPath, sep = '\t', comment = '#', names = colnames)
samples = VCF.columns.tolist()[9:VCF.shape[0]]
# Add Type
VCF['TYPE'] = [x.split('_')[1] for x in VCF['ID']]
# Filter for Type if required
if types:
	logging.info(f'Filter VCF for {" ".join(types)}')
	VCF = VCF.loc[VCF['TYPE'].isin(types),:]
# Add Len
VCF['LEN'] = [abs(int(x.split('_')[5])) for x in VCF['ID']]

# Add second breakpoint for each SV:
## Insertions: a single breakpoint (BP1 == BP2)
VCF.loc[VCF['TYPE'] == 'INS', 'POS2'] = VCF.loc[VCF['TYPE'] == 'INS', 'POS']
VCF.loc[VCF['TYPE'] == 'INS', 'CHROM2'] = VCF.loc[VCF['TYPE'] == 'INS', '#CHROM']
## Deletion/Duplication/Contraction/Inversion counts for Pos and Pos+Len
VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'POS2'] = VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'POS'] + VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'LEN']
VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'CHROM2'] = VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), '#CHROM']
## Translocations count for 2 breakpoints
VCF.loc[VCF['TYPE'] == 'TRA', 'POS2'] = [int(x.split(';')[-1].split(':')[-1]) for x in VCF.loc[VCF['TYPE'] == 'TRA', 'INFO']]
VCF.loc[VCF['TYPE'] == 'TRA', 'CHROM2'] = [x.split(';')[-1].split(':')[0].split('=')[1] for x in VCF.loc[VCF['TYPE'] == 'TRA', 'INFO']]
## Rename first breakpoint
VCF = VCF.rename(columns = {'#CHROM': 'CHROM1', 'POS': 'POS1'})

# Split by haplotypes
haploVCF = []
for s in samples:
	haploVCF += [VCF[s].str.split('|', expand = True).rename(columns = {0:f'{s}1', 1:f'{s}2'})]
haploVCF = pd.concat(haploVCF, axis = 1)
haplotypes = haploVCF.columns.tolist()

# Extract windows over which Pi and Theta must be calculated ==================
logging.info('Windows definition')

windows = []
# Extract chromosomes in VCF and length
chromosomes = [(x.split(',')[0].split('=')[2], int(x.split('=')[3].split('>')[0])) for x in header if x.startswith('##contig=<ID=')]
for c, size in chromosomes:
	# iterate on each window
	max_start = size + sliding - window
	# When window is larger than chromosome size, set max_start to 1 (to have a single wondow on the chromosome)
	if max_start < 0:
		max_start = 1
	for s in range(0, max_start, sliding): 
		e = min(s + window, size)
		windows += [pd.DataFrame(data = {'Chromosome': [c], 'Start': [s], 'End': [e]})]
windows = pd.concat(windows).reset_index(drop = True)

# Pi Definition ===============================================================
def computePi(index):
	'''Compute Structural Theta Pi over a window defined by its index in teh windows data frame'''
	# Get location of window
	c = windows.loc[index, 'Chromosome']
	s = windows.loc[index, 'Start']
	e = windows.loc[index, 'End']
	# Filter SVs in window and join haplotypes
	subHaplo = haploVCF.loc[
		((VCF['CHROM1'] == c) & (VCF['POS1'].isin(list(range(s, e+1))))) | # First BP is in window
		((VCF['CHROM2'] == c) & (VCF['POS2'].isin(list(range(s, e+1))))) | # Second BP is in window
		((VCF['TYPE'] != 'TRA') & (VCF['CHROM1'] == c) & (VCF['CHROM2'] == c) & (VCF['POS1'] <= s) & (VCF['POS2'] >= e)) # Window is between both BP
		,:].apply(''.join)
	# Count presence of each haplotype
	subHaplo = subHaplo.value_counts()
	subHaplo = pd.DataFrame(data = {'Haplotype': subHaplo.index, 'Count': subHaplo}).reset_index(drop = True)
	subHaplo['Freq'] = subHaplo['Count'] / len(haplotypes)
	# Generate all possible combination of haplotypes
	comb = list(itertools.combinations(subHaplo.index.tolist(), 2))
	comb_elem1 = [x[0] for x in comb]
	comb_elem2 = [x[1] for x in comb]
	combDF = pd.DataFrame(data = {'HP1': subHaplo.loc[comb_elem1, 'Haplotype'].reset_index(drop = True), 
								'HP2': subHaplo.loc[comb_elem2, 'Haplotype'].reset_index(drop = True), 
								'Freq1': subHaplo.loc[comb_elem1, 'Freq'].reset_index(drop = True), 
								'Freq2': subHaplo.loc[comb_elem2, 'Freq'].reset_index(drop = True)})
	# Count number of differences for each combination
	# The number of differences is divided by two to compensate for expansion of the VCF above
	combDF['NDiff'] = [sum(1 for a, b in zip(combDF.loc[i,'HP1'], combDF.loc[i,'HP2']) if a != b) for i in combDF.index]
	# Compute Pi (Number of differences / length of window)
	# Compute Factor that will be summed to obtain average Structural Pi on the window (2*xi*xj*πij)
	combDF['Pi'] = combDF['NDiff'] / window
	combDF['Factor'] = 2 * combDF['Freq1'] * combDF['Freq2'] * combDF['Pi']
	# Compute Structural diversity on the window
	Pi = (len(haplotypes) / (len(haplotypes) - 1)) * combDF['Factor'].sum()
	# Return DataFrame with Result
	return(pd.DataFrame(data = {'Chromosome': [c], 'Start': [s], 'End': [e], 'Statistics': ['ThetaPi'], 'Value': [Pi]}))

# ThetaS definition ============================================================
def computeTheta(index):
	'''Compute Structural ThetaS over a window defined by its index in the windows data frame'''
	# Get location of window
	c = windows.loc[index, 'Chromosome']
	s = windows.loc[index, 'Start']
	e = windows.loc[index, 'End']
	# Filter SVs in window
	subHaplo = haploVCF.loc[
		((VCF['CHROM1'] == c) & (VCF['POS1'].isin(list(range(s, e+1))))) | # First BP is in window
		((VCF['CHROM2'] == c) & (VCF['POS2'].isin(list(range(s, e+1))))) | # Second BP is in window
		((VCF['TYPE'] != 'TRA') & (VCF['CHROM1'] == c) & (VCF['CHROM2'] == c) & (VCF['POS1'] <= s) & (VCF['POS2'] >= e)) # Window is between both BP
		,:]
	# Count number of SVs in the window
	NbSVs = subHaplo.shape[0]
	# Compute ThetaS
	ThetaS = NbSVs / harmonicNumber(len(haplotypes) - 1) / window
	# Return DataFrame with Result
	return(pd.DataFrame(data = {'Chromosome': [c], 'Start': [s], 'End': [e], 'Statistics': ['ThetaS'], 'Value': [ThetaS]}))

# Gene diversity Definition ===============================================================
def computeGeneDiv(index):
	'''Compute Structural Gene diversity (See Harris et al., 2017, G3) over a window defined by its index in the windows data frame'''
	# Get location of window
	c = windows.loc[index, 'Chromosome']
	s = windows.loc[index, 'Start']
	e = windows.loc[index, 'End']
	# Filter SVs in window and join haplotypes
	subHaplo = haploVCF.loc[
		((VCF['CHROM1'] == c) & (VCF['POS1'].isin(list(range(s, e+1))))) | # First BP is in window
		((VCF['CHROM2'] == c) & (VCF['POS2'].isin(list(range(s, e+1))))) | # Second BP is in window
		((VCF['TYPE'] != 'TRA') & (VCF['CHROM1'] == c) & (VCF['CHROM2'] == c) & (VCF['POS1'] <= s) & (VCF['POS2'] >= e)) # Window is between both BP
		,:].apply(''.join)
	# Count presence of each haplotype
	subHaplo = subHaplo.value_counts()
	subHaplo = pd.DataFrame(data = {'Haplotype': subHaplo.index, 'Count': subHaplo}).reset_index(drop = True)
	subHaplo['Freq'] = subHaplo['Count'] / len(haplotypes)
	subHaplo['FreqSquare'] = [x ** 2 for x in subHaplo['Freq']]
	# Compute Gene diversity
	Het = (len(haplotypes) / (len(haplotypes) - 1)) * ( 1 - subHaplo['FreqSquare'].sum() )
	# Return DataFrame with Result
	return(pd.DataFrame(data = {'Chromosome': [c], 'Start': [s], 'End': [e], 'Statistics': ['GeneDiversity'], 'Value': [Het]}))


# Run computation in parallel =================================================

# Create the pool
process_pool = Pool(processes=nthreads)
# Start processes in the pool
logging.info(f'Compute Pi over {windows.shape[0]} windows')
outDF = process_pool.map(computePi, windows.index.tolist())
logging.info(f'Compute Theta over {windows.shape[0]} windows')
outDF += process_pool.map(computeTheta, windows.index.tolist())
logging.info(f'Compute Gene Diversity over {windows.shape[0]} windows')
outDF += process_pool.map(computeGeneDiv, windows.index.tolist())

outDF = pd.concat(outDF).reset_index(drop = True)
outDF.to_csv(outputPath, sep = '\t', index = False)



