#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/10/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script changes the ID of the variants in a VCF file. It requires a
tab-separated file with col1: old ID and col2: new ID. 
'''
# ---------------------------------------------------------------------------
import gzip
import pandas as pd
import argparse
import os
import logging
import sys
sys.path.insert(1, '/home/vloegler/SVCalling_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
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
parser.add_argument("-id", "--id", help="Tab-separated file with old IDs (col1) and new IDs (col2). No header. ", required=True)
parser.add_argument("-o", "--output", help="Output file (VCF format, not gzipped)", required=True)

# Read arguments from the command line
args = parser.parse_args()

vcfPath=os.path.abspath(args.vcf)
idPath=os.path.abspath(args.id)
outputPath=os.path.abspath(args.output)

#vcfPath = 'SV.1087Samples.vcf'
#vcfPath = 'test.vcf.gz'
#idPath='IDs.txt'

# Read VCF =============================================================

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
VCF.index = VCF['ID'].copy()

# Change IDs ===========================================================

# Read ID table
IDs = pd.read_csv(idPath, sep = '\t', names = ['OldID', 'NewID'])
IDs.index = IDs['OldID']

# Remove IDs absent from the VCF
IDs = IDs.loc[IDs['OldID'].isin(VCF['ID'])]

# Change IDs
VCF.loc[IDs['OldID'],'ID'] = IDs['NewID']
logging.info(f'{IDs.shape[0]} sequence IDs changed. ')

# Output file ==========================================================
VCF = VCF.rename(columns={"ALT_ALLELE": "ALT"})

# Write header
with open(outputPath, 'w') as out:
        for line in header:
            out.write(line)

VCF.to_csv(outputPath, sep = '\t', mode='a', index = False)
