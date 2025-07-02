#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/10/10
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script find SV-gene pairs using a SV vcf and a GFF file. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import argparse
import os
import sys
sys.path.insert(1, '/home/vloegler/SVCalling_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
import gzip
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
parser.add_argument("-v", "--vcf", help="Path of the SV VCF (vcf or vcf.gz)", required=True)
parser.add_argument("-g", "--gff", help="Path to the annotation file (GFF3)", required=True)
parser.add_argument("-k", "--kb", help="Size (kb) of the regions considered as impacted around SV breakpoints", required=True, type = int)
parser.add_argument("-o", "--output", help="Output file (VCF format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

vcfPath=os.path.abspath(args.vcf)
gffPath=os.path.abspath(args.gff)
margin = args.kb
outputPath=os.path.abspath(args.output)

#vcfPath = '04_PopulationVCF/SV.1087Samples.vcf.gz'
#gffPath='../03-data/sace_R64-3-1_20210421.gff'

# Get region impacted by each SV ==============================================
logging.info('Get region impacted by each SV')

# Get Column names of VCF
if vcfPath.endswith('.gz'):
    with gzip.open(vcfPath, "rb") as vcf:
        for line in vcf:
            line = line.decode()
            if line.startswith('#CHROM'):
                colnames = line.strip().split('\t')
                break
else:
    with open(vcfPath, "r") as vcf:
        for line in vcf:
            if line.startswith('#CHROM'):
                colnames = line.strip().split('\t')
                break

# Change ALT colnames to ALT ALLELE to avoid confonding with ALT strain
colnames[4] = 'ALT_ALLELE'

# Read vcf
VCF = pd.read_csv(vcfPath, sep = '\t', comment = "#", names = colnames)
VCF = VCF.loc[:,['ID', '#CHROM', 'POS', 'INFO']]
VCF.loc[:,'TYPE'] = [x.split('_')[1] for x in VCF['ID']]

# Split VCF by type and get breakpoint: INS (1 breakpoint), DEL-DUP-CONTR-INV (2 breakpoints on the same chromosome), TRA (2 breakpoints)
## Insertions
VCF_INS = VCF.loc[VCF['TYPE'] == 'INS', ['ID', 'TYPE', '#CHROM', 'POS']].copy()
VCF_INS = VCF_INS.rename(columns = {'#CHROM': 'CHROM'})
## Deletions Duplications Contractions Inversions (DDCI)
VCF_DDCI = VCF.loc[VCF['TYPE'].isin(['DEL', 'DUP', 'CONTR', 'INV']), ['ID', 'TYPE', '#CHROM', 'POS', 'INFO']].copy()
VCF_DDCI = VCF_DDCI.rename(columns = {'#CHROM': 'CHROM1', 'POS': 'POS1'})
VCF_DDCI.loc[:,'LEN'] = [abs(int([x for x in INFO.split(';') if x.startswith('SVLEN')][0].split('=')[1])) for INFO in VCF_DDCI['INFO']]
VCF_DDCI.loc[:,'CHROM2'] = VCF_DDCI['CHROM1']
VCF_DDCI.loc[:,'POS2'] = VCF_DDCI['POS1'] + VCF_DDCI['LEN']
VCF_DDCI = VCF_DDCI.loc[:,['ID', 'TYPE', 'CHROM1', 'POS1', 'CHROM2', 'POS2']]
## Translocations
VCF_TRA = VCF.loc[VCF['TYPE'] == 'TRA', ['ID', 'TYPE', '#CHROM', 'POS', 'INFO']].copy()
VCF_TRA = VCF_TRA.rename(columns = {'#CHROM': 'CHROM1', 'POS': 'POS1'})
VCF_TRA.loc[:,'BP2'] = [[x for x in INFO.split(';') if x.startswith('BREAKPOINT2')][0].split('=')[1] for INFO in VCF_TRA['INFO']]
VCF_TRA.loc[:,'CHROM2'] = [x.split(':')[0] for x in VCF_TRA['BP2']]
VCF_TRA.loc[:,'POS2'] = [int(x.split(':')[1]) for x in VCF_TRA['BP2']]
VCF_TRA = VCF_TRA.loc[:,['ID', 'TYPE', 'CHROM1', 'POS1', 'CHROM2', 'POS2']]

# Add bed object of SVs (with X kb margins around SVs)

## Insertions: X kb upstream and downstream breakpoint
VCF_INS.loc[:,'BED'] = [BED(BEDcoordinates(VCF_INS.loc[i,'CHROM'], max(0, VCF_INS.loc[i,'POS'] - 1 - margin*1000), VCF_INS.loc[i,'POS'] + margin*1000)) for i in VCF_INS.index]
VCF_INS = VCF_INS.loc[:,['ID', 'TYPE', 'BED']]

## Deletions Duplications Contractions: X kb upstream and downstream breakpoints, plus region in SV
VCF_DDC = VCF_DDCI.loc[VCF_DDCI['TYPE'].isin(['DEL', 'DUP', 'CONTR']),:].copy()
VCF_DDC.loc[:,'BED'] = [BED(BEDcoordinates(VCF_DDC.loc[i,'CHROM1'], max(0, VCF_DDC.loc[i,'POS1'] - 1 - margin*1000), VCF_DDC.loc[i,'POS2'] + margin*1000)) for i in VCF_DDC.index]
VCF_DDC = VCF_DDC.loc[:,['ID', 'TYPE', 'BED']]

## Inversions and translocations: X kb upstream and downstream breakpoints, withous region in SV
VCF_INV = VCF_DDCI.loc[VCF_DDCI['TYPE'] == 'INV'].copy()
VCF_IT = pd.concat([VCF_TRA, VCF_INV])
VCF_IT.loc[:,'BED1'] = [BED(BEDcoordinates(VCF_IT.loc[i,'CHROM1'], max(0, VCF_IT.loc[i,'POS1'] - 1 - margin*1000), VCF_IT.loc[i,'POS1'] + margin*1000)) for i in VCF_IT.index]
VCF_IT.loc[:,'BED2'] = [BED(BEDcoordinates(VCF_IT.loc[i,'CHROM2'], max(0, VCF_IT.loc[i,'POS2'] - 1 - margin*1000), VCF_IT.loc[i,'POS2'] + margin*1000)) for i in VCF_IT.index]
VCF_IT.loc[:,'BED'] = [VCF_IT.loc[i,'BED1'] + VCF_IT.loc[i,'BED2'] for i in VCF_IT.index]
VCF_IT = VCF_IT.loc[:,['ID', 'TYPE', 'BED']]

# Merge everything
SVdata = pd.concat([VCF_INS, VCF_DDC, VCF_IT])


# Read GFF file ===============================================================
gff = pd.read_csv(gffPath, sep = '\t', comment = '#', names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
gff = gff.loc[(gff['type'] == 'gene') & (gff['seqid'] != 'chrmt'),:]
gff.loc[:,'GeneName'] = [[x.split('=')[1] for x in att.split(';') if x.startswith('Name=')][0] for att in gff['attributes']]
gff.loc[:,'Gene'] = [[x.split('=')[1] for x in att.split(';') if x.startswith('gene=')][0] if 'gene=' in att else 'NA' for att in gff['attributes']]
gff.loc[:,'BED'] = [BED(BEDcoordinates(gff.loc[i,'seqid'], gff.loc[i,'start'] - 1, gff.loc[i,'end'])) for i in gff.index]
gff = gff.reset_index(drop = True)
gff['ID'] = gff.index

# Compute overlap between each combination of SV and Genes ====================

## Add one column for each gene with the fraction of the gene covered by the SV
logging.info('Compute overlap between SVs and genes')
overlapSeries = []
SVdata['Chrs'] = [x.IDs for x in SVdata['BED']]
for i in gff.index:
	geneChr = gff.loc[i, 'seqid']
	subData = SVdata.loc[[geneChr in x for x in SVdata['Chrs']],:].copy() # Do not check overlap for SVs on a different chromosome
	subData['BED2'] = gff.loc[i, 'BED']
	subData['GeneLen'] = len(gff.loc[i, 'BED'])
	subData['NonAffectedPart'] = subData['BED2'] - subData['BED']
	subData['NonAffectedPart'] = subData['NonAffectedPart'].apply(len)
	subData['AffectedPart'] = subData['GeneLen'] - subData['NonAffectedPart']
	GeneCol = subData['AffectedPart'] / subData['GeneLen']
	GeneCol.name = f'Gene{i}'
	overlapSeries += [GeneCol]

SVdata = SVdata.drop(columns = ['Chrs'])
data = pd.concat([SVdata] + overlapSeries, axis = 1)

logging.info('Process results')
# Pivot longer and drop SV gene pairs with no overlap
data = pd.wide_to_long(data, stubnames = 'Gene', i = ['ID', 'TYPE', 'BED'], j = 'GeneID')
data = data.reset_index()
data = data.rename(columns = {'Gene': 'Overlap'})
data = data.loc[data['Overlap'] > 0,:]
data['GeneName'] = [gff.loc[i,'GeneName'] for i in data['GeneID']]
data['Gene'] = [gff.loc[i,'Gene'] for i in data['GeneID']]

# Order by position on genome
data['N_SV'] = [int(id.split('_')[0]) for id in data['ID']]
data = data.sort_values(by = 'N_SV')

# Drop unwanted columns
data = data.drop(columns = ['TYPE', 'BED', 'GeneID', 'N_SV'])
data = data.loc[:,['ID', 'GeneName', 'Gene', 'Overlap']]
data = data.rename(columns = {'ID': 'SV_ID'})

# Output table ================================================================
data.to_csv(outputPath, sep = '\t', index = False)

