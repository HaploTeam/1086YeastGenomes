#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/10/04
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script filter SV haplotypes for a defined region of the genome. 
'''
# ---------------------------------------------------------------------------
import gzip
import pandas as pd
import argparse
import os
import logging
logging.basicConfig(level=logging.NOTSET)
# ---------------------------------------------------------------------------
# Definitions

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="Path of the VCF file (gzipped or not)", required=True)
parser.add_argument("-c", "--chromosome", help="Chromosome ID", required=True)
parser.add_argument("-s", "--start", help="Start position of the region", required=True, type = int)
parser.add_argument("-e", "--end", help="End position of the region", required=True, type = int)
parser.add_argument("-o", "--output", help="Output file (tsv format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

vcfPath=os.path.abspath(args.vcf)
chromosome = args.chromosome
start = args.start
end = args.end
outputPath=os.path.abspath(args.output)

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
# Add Len
VCF['LEN'] = [abs(int(x.split('_')[5])) for x in VCF['ID']]

# Add second breakpoint information
## Insertions have same than first bp
VCF.loc[VCF['TYPE'] == 'INS', 'POS2'] = VCF.loc[VCF['TYPE'] == 'INS', 'POS']
VCF.loc[VCF['TYPE'] == 'INS', 'CHROM2'] = VCF.loc[VCF['TYPE'] == 'INS', '#CHROM']
## Deletion/Duplication/Contraction/Inversion counts for Pos and Pos+Len
VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'POS2'] = VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'POS'] + VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'LEN']
VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), 'CHROM2'] = VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), '#CHROM']
## Translocations count for 2 breakpoints
VCF.loc[VCF['TYPE'] == 'TRA', 'POS2'] = [int(x.split(';')[-1].split(':')[-1]) for x in VCF.loc[VCF['TYPE'] == 'TRA', 'INFO']]
VCF.loc[VCF['TYPE'] == 'TRA', 'CHROM2'] = [x.split(';')[-1].split(':')[0].split('=')[1] for x in VCF.loc[VCF['TYPE'] == 'TRA', 'INFO']]
## Rename columns
VCF = VCF.rename(columns = {'#CHROM': 'CHROM1', 'POS': 'POS1'})

# Split by haplotypes
haploVCF = []
for s in samples:
	haploVCF += [VCF[s].str.split('|', expand = True).rename(columns = {0:f'{s}1', 1:f'{s}2'})]
haploVCF += [VCF['ID']]
haploVCF = pd.concat(haploVCF, axis = 1)
haplotypes = haploVCF.drop(columns = ['ID']).columns.tolist()


# Filter region =====================================================

# Extract SVs located in the window
subHaploVCF = haploVCF.loc[
	((VCF['CHROM1'] == chromosome) | (VCF['CHROM2'] == chromosome)) & # Check for chromosome
	(VCF['POS1'].isin(list(range(start, end+1))) | # First BP is in window
		VCF['POS2'].isin(list(range(start, end+1))) | # Second BP is in window
		((VCF['CHROM1'] == VCF['CHROM2']) & (VCF['POS1'] <= start) & (VCF['POS2'] >= end))) # Window is between both BP
	,:]

# Get haplotypes and frequencies
out = subHaploVCF.drop(columns = ['ID']).apply(';'.join).value_counts()
out = pd.DataFrame(out).reset_index().rename(columns = {'index': 'Haplotype', 'count': 'HaplotypeCount'})
out[subHaploVCF['ID']] = out['Haplotype'].str.split(';',expand=True)
out = out.copy()
out = out.drop(columns = ['Haplotype'])
out['Chromosome'] = chromosome
out['Start'] = start
out['End'] = end
out = out.loc[:,['Chromosome', 'Start', 'End', 'HaplotypeCount'] + subHaploVCF['ID'].tolist()]

# Output file
out.to_csv(outputPath, sep = '\t', index = False)
