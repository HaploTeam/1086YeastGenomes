#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/06/20
# version ='1.1'
# ---------------------------------------------------------------------------
'''
This script merges all the information on the pangenome non-reference genes 
in a single table. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import argparse
import os
import sys
sys.path.insert(1, '/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
# ---------------------------------------------------------------------------
# Definitions

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-l", "--locationTable", help="Localization table (output from GeneClustering.py)", required=True)
parser.add_argument("-p", "--nonRedundantPangenome", help="Fasta file of the pangenome without redundancy (output from RemoveRedundancy.py)", required=True)
parser.add_argument("-r", "--originTable", help="Origin table (output from InferGeneOrigin.py)", required=True)
parser.add_argument("-f", "--functionTable", help="Function table (output from InferGeneFunction.py)", required=True)
parser.add_argument("-i", "--interproscanTable", help="tsv output of interproscan", required=True)
parser.add_argument("-o", "--output", help="Output tsv file", required=True)

# Read arguments from the command line
args = parser.parse_args()

locationTable = os.path.abspath(args.locationTable)
NonRedundantFasta = os.path.abspath(args.nonRedundantPangenome)
originTable = os.path.abspath(args.originTable)
functionTable = os.path.abspath(args.functionTable)
interproscanTable = os.path.abspath(args.interproscanTable)
output = os.path.abspath(args.output)

locationTable = '/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains/Pangenome.NonRef.Localization.tsv'
originTable = '/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/09_OriginInference/Pangenome.origin.tsv'
functionTable = '/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/10_FunctionInference/Pangenome.FunctionTransfer.tsv'
interproscanTable = '/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/11_InterProScan/Pangenome.NoRedundancy.pep.faa.InterProScan.tsv'
NonRedundantFasta = '/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.cds.fna'

# Import tables

## Location
locationTable = pd.read_csv(locationTable, sep = '\t')
locationTable.index = locationTable['GeneID']

## Phylogenetic origin
originTable = pd.read_csv(originTable, sep = '\t')
originTable.index = originTable['qseqid']

## Function
functionTable = pd.read_csv(functionTable, sep = '\t')
functionTable.index = functionTable['qseqid']

## Interproscan
interproscanTable = pd.read_csv(interproscanTable, sep = '\t', names = ['GeneID', 'MD5', 'Length', 'Analysis', 'SignatureAccession', 'SignatureDescription', 'Start', 'Stop', 'Score', 'Status', 'Date', 'InterProAccession', 'InterProDescription', 'GOannotation', 'PathwaysAnnotations'])
interproscanTable = interproscanTable.loc[:,['GeneID', 'GOannotation']]
interproscanTable = interproscanTable.loc[interproscanTable['GeneID'].isin(locationTable['GeneID']),:]
interproscanTable['GOannotation'] = [[y.split('(')[0] for y in x.split('|')] if x != '-' else [] for x in interproscanTable['GOannotation']]
interproscanTable = interproscanTable.groupby('GeneID').agg({'GOannotation':'sum'})
interproscanTable['GOannotation'] = [' '.join(set(sorted(x))) for x in interproscanTable['GOannotation']]

# Add redundancy information
locationTable['RedundantGene'] = True
locationTable.loc[locationTable['GeneID'].isin(Fasta(NonRedundantFasta).getID()), 'RedundantGene'] = False
locationTable = locationTable.loc[:,['GeneID', 'RedundantGene', 'RepresentantGeneID', 'OriginStrain', 'Chromosome', 'Strand', 'Start', 'End', 'PreviousGene', 'NextGene']]

# Merge tables
final = pd.concat([locationTable, functionTable.loc[:,['Final', 'TransferredGOterms']], interproscanTable, originTable.loc[:,['origin', 'mechanism', 'sscinames', 'sseqid']]], axis = 1)
final = final.rename(columns={"Final": "Function", "GOannotation":"InterProScan_GO_Terms", "origin": "Origin", "mechanism":"Mechanism", "sscinames": "Species", "sseqid": "OriginProtein"})
final.loc[final['TransferredGOterms'].isna(), 'TransferredGOterms'] = ''
final.loc[final['InterProScan_GO_Terms'].isna(), 'InterProScan_GO_Terms'] = ''
final.loc[:,'Total_GO_Terms'] = [' '.join(filter(lambda x: x != '', set(sorted(final.loc[i,'TransferredGOterms'].split(' ') + final.loc[i,'InterProScan_GO_Terms'].split(' '))))) for i in final.index]

# Write table
final.to_csv(output, sep = '\t', index = False)
