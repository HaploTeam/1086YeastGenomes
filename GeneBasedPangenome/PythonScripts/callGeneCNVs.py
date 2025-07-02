#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/09/08
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script use Bam mapping file of short reads on gene-based pangenome to 
estimate the number of copy of each gene of the pangenome in an isolate. 
'''
# ---------------------------------------------------------------------------
import os
import argparse
import subprocess
import pandas as pd
import numpy as np
# ---------------------------------------------------------------------------
# Definitions

# ---------------------------------------------------------------------------

# =============
# GET ARGUMENTS
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = '''
This script use Bam mapping file of short reads on gene-based pangenome to 
estimate the number of copt of each gene of the pangenome in an isolate. 
'''
)
#parser.add_argument("-f", "--fastaPangenome", help="Fasta file of the pangenome", required=True)
parser.add_argument("-B", "--BAM", help="Mapping file of short reads sequencing on the pangenome (BAM format)", required=True)
parser.add_argument("-s", "--strain", help="Strain name", required=True)
parser.add_argument("-p", "--ploidy", help="Ploidy of the strain", type = int, required=True)
parser.add_argument("-i", "--identicalRegions", help="Tsv file of the identical regions in the pangenome", required=True)
parser.add_argument("-x", "--x", help="x value to be used for gene presence detection (Presence is proportional to NormDepth * Ploidy^x)", type = float, required = True)
parser.add_argument("-t", "--threshold", help="Minimum value of (NormDepth * Ploidy^x) required for detection of gene presence", type = float, default = 0.5)
parser.add_argument("-o", "--output", help="Name of the output file", required=True)

# Read arguments from the command line
args = parser.parse_args()
#fastaPath = os.path.abspath(args.fastaPangenome)
bamPath = os.path.abspath(args.BAM)
strain = args.strain
ploidy = int(args.ploidy)
identicalRegionPath = os.path.abspath(args.identicalRegions)
x = args.x
threshold = args.threshold
outputPath = os.path.abspath(args.output)

#bamPath = "1_Mapping/AAA.bam"
#identicalRegionPath = "../03-data/Pangenome.NoRedundancy.cds.IdenticalRegions.tsv"
#ploidy = 1

# Get coverage depth along pangenome
command = ['samtools', 'depth', '-aa', bamPath, '-o', f'{bamPath}.depth']
subprocess.call(command)
covDF = pd.read_csv(f'{bamPath}.depth', sep = "\t", names = ['Gene', 'Pos', 'Depth'])
covDFwide = covDF.pivot(index = "Gene", columns = "Pos")
covDFwide.columns = covDFwide.columns.droplevel(0)

# Read identical regions bedfile
idenDF = pd.read_csv(identicalRegionPath, sep = "\t")
# Mask identical regions coverage
for i in idenDF.index:
    covDFwide.loc[idenDF.loc[i,'Gene1'], list(range(idenDF.loc[i,'Start1']+1, idenDF.loc[i,'End1']+1))] = np.nan
# Pivot DF longer
covDF = pd.melt(covDFwide, ignore_index = False).reset_index()
covDF.columns = ['Gene', 'Pos', 'Depth']
covDF = covDF.dropna(how = 'any')


# Compute normalized coverage for each gene

### Get median coverage on the whole pangenome
medCov = covDF.loc[(covDF['Depth'] > 0), 'Depth'].median() # Median coverage of the pangenome
### Get median coverage of each gene
covDF = covDF.groupby('Gene').agg({'Depth': 'median'}).reset_index()
### Normalize by median coverage along the pangenome
covDF.loc[:,'NormDepth'] = covDF['Depth'] / medCov
### Add ploidy
covDF['Ploidy'] = ploidy
## Add the NormDepth * Ploidy^x column to detect gene presence
covDF['NormDepthTimesPloidyX'] = [covDF.loc[i,'NormDepth'] * (covDF.loc[i,'Ploidy']**x) for i in covDF.index]
### Add the NormDepth * Ploidy column to estimate the number of copies present
covDF['NormDepthTimesPloidy'] = covDF['NormDepth'] * covDF['Ploidy']

# Determine presence and number of copies for each gene

### Determine presence of gene
covDF['Presence'] = covDF['NormDepthTimesPloidyX'] >= threshold
### Determine number of copies per haplotype
covDF['CN'] = round(covDF['NormDepthTimesPloidy']) / ploidy
### Set number of copies to 1 when gene was detected as present and CN is 0
covDF.loc[(covDF['Presence'] == True) & (covDF['CN'] == 0), 'CN'] = 1 / covDF.loc[(covDF['Presence'] == True) & (covDF['CN'] == 0), 'Ploidy']
### Set number of copies to 0 when gene was detected as absent
covDF.loc[covDF['Presence'] == False, 'CN'] = 0

# Output file
covDF['Strain'] = strain
covDF = covDF.loc[:,['Strain', 'Ploidy', 'Gene', 'Depth', 'NormDepth', 'CN']]
covDF.to_csv(outputPath, sep = '\t', index = False)

