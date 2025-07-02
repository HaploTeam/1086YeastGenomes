#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2024/02/28
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script takes a plink file and output a list of SNPs, INDELs and SVs by 
looking at the alleles in the bim file. 
'''
# ---------------------------------------------------------------------------
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
parser = argparse.ArgumentParser('''
This script takes a plink file and output a list of SNPs, INDELs and SVs by 
looking at the alleles in the bim file. 
''')
parser.add_argument("-b", "--bfile", help="Plink file", required=True)

# Read arguments from the command line
args = parser.parse_args()

bfile=os.path.abspath(args.bfile)

bfile = '/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis/11_HeritabilityEstimates/GraphGenotyping.937Samples.Reheader.DP2.TrimAlt.Atomize.VariantsAnnotated.plink'

# Read bim file
bim = pd.read_csv(f'{bfile}.bim', sep = '\t| ', names = ['Chrom', 'ID', 'Dist', 'Pos', 'Major', 'Minor'], engine='python')
# Compute length of alleles
bim['MajorLen'] = bim['Major'].apply(len)
bim['MinorLen'] = bim['Minor'].apply(len)
# Deduce type of allele
bim['Type'] = 'SV'
bim.loc[(bim['MajorLen'] < 50) & (bim['MinorLen'] < 50), 'Type'] = 'INDEL'
bim.loc[bim['MajorLen'] == bim['MinorLen'], 'Type'] = 'SNP'

# Write list of SNPs, INDELs and SVs
bim.loc[bim['Type'] == 'SNP', 'ID'].to_csv(f'{bfile}.snps', header = False, index = False)
bim.loc[bim['Type'] == 'INDEL', 'ID'].to_csv(f'{bfile}.indels', header = False, index = False)
bim.loc[bim['Type'] == 'SV', 'ID'].to_csv(f'{bfile}.svs', header = False, index = False)
