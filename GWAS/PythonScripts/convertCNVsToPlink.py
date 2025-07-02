#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2024/01/03
# version ='1.1'
# ---------------------------------------------------------------------------
'''
This script converts a tsv files with CNV data (with columns Strain, Ploidy, 
Gene, Depth, NormDepth, CN) into a plink ped and map files. The NormDetph 
column correponds to the coverage of the median coverage of the gene over the 
median coverage of the pangenome. 
Each variant correspond to a gene and a NormDepth threshold. A value of 1 
means that the NormCoverage is above the threshold, a value of 2 is the 
opposite. 
Ex:
Normdepth = 0.42
Variants    GenotypeHP1 GenotypeHP2
Gene1_Cov0.25      1           1
Gene1_Cov0.75      2           2

Normdepth = 1.02
Variants    GenotypeHP1 GenotypeHP2
Gene1_Cov0.25      1           1
Gene1_Cov0.75      1           1
Gene1_Cov1.25      2           2
'''
# ---------------------------------------------------------------------------
import os
import re
import csv
import argparse
import pandas as pd
import math
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
This script converts a tsv files with CNV data (with columns Strain, Ploidy, 
Gene, Depth, NormDepth, CN) into a plink ped and map files. The NormDetph 
column correponds to the coverage of the median coverage of the gene over the 
median coverage of the pangenome. 
Each variant correspond to a gene and a NormDepth threshold. A value of 1 
means that the NormCoverage is above the threshold, a value of 2 is the 
opposite. 
Ex:
Normdepth = 0.42
Variants    GenotypeHP1 GenotypeHP2
Gene1_Cov0.25      1           1
Gene1_Cov0.75      2           2

Normdepth = 1.02
Variants    GenotypeHP1 GenotypeHP2
Gene1_Cov0.25      1           1
Gene1_Cov0.75      1           1
Gene1_Cov1.25      2           2
'''
)

parser.add_argument("-i", "--input", help="CNV Tsv file with header. Tsv file must at least contain 3 columns: Strain, Gene and CN (which refers to the gene copy number normalized by ploidy). ", required=True)
parser.add_argument("-g", "--gff", help="GFF file containing gene location information", required=True)
parser.add_argument("-l", "--additionalLocation", help="Additional location file for non reference genes (Pangenome.Annotation.tsv file). Columns should be :GeneID, RedundantGene, RepresentantGeneID, OriginStrain, Chromosome, Strand, Start, End, PreviousGene, NextGene, Function, TransferredGOterms, InterProScan_GO_Terms, Origin, Mechanism, Species, OriginProtein, ProteInfer_GO_Terms, ProteInfer_PFAM_Domains, ProteInfer_EC_Numbers, Total_GO_Terms; although only GeneID, Chromosome and Start will be considered. ", required=True)
parser.add_argument("-o", "--output", help="Prefix of the output plink matrix (ped and map files)", required=True)

# Read arguments from the command line
args = parser.parse_args()
inputPath = os.path.abspath(args.input)
gffPath = os.path.abspath(args.gff)
additionalLocationPath = os.path.abspath(args.additionalLocation)
outputPath = os.path.abspath(args.output)

#inputPath='/ccc/scratch/cont007/fg0006/loeglerv/GWAS_CNV/Genotypes/full1086Sace_MatrixCNV.tsv.gz'
#gffPath='/ccc/scratch/cont007/fg0006/loeglerv/GWAS_CNV/Genotypes/Sace_SGD_R64-4-1_20230830.gff'
#additionalLocationPath='/ccc/scratch/cont007/fg0006/loeglerv/GWAS_CNV/Genotypes/Pangenome.NonRef.Localization.tsv'

# Import data
# ===========
# Get location information from gff and additional location file
gff = pd.read_csv(gffPath, sep = '\t', comment = '#', names = ['seqid', 'source', 'type', 'start', 'end', 'score', 'strand', 'phase', 'attributes'])
gff = gff.loc[gff['type'].isin(['gene', 'pseudogene', 'transposable_element_gene', 'blocked_reading_frame']),:]
gff['Gene'] = gff['attributes'].apply(lambda x: x.split(';')[0].split('=')[1])
gff = gff.loc[gff['seqid'] != 'chrmt', ['Gene', 'seqid', 'start']].rename(columns = {'seqid': 'Chromosome', 'start': 'Start'})

addLoc = pd.read_csv(additionalLocationPath, sep = '\t', header = 0)
addLoc = addLoc.loc[:,['GeneID', 'Chromosome', 'Start']].rename(columns = {'GeneID': 'Gene'})
def getChromosome(x):
    if x.count('chromosome') == 1:
        return([y.split('.')[0] for y in x.split('_') if 'chromosome' in y][0])
    else:
        return('unplaced')

addLoc['Chromosome'] = addLoc['Chromosome'].apply(getChromosome)
location = pd.concat([gff, addLoc], axis = 0)
# Set chromosome as numeric
# Unplaced genes are set to chromosome 17
location['Chromosome'] = [int(re.sub('chromosome', '', x)) if 'chromosome' in x else 17 for x in location['Chromosome']]

# Import Coverage information
data = pd.read_csv(inputPath, sep = '\t', header = 0)
data = data.loc[:,['Strain', 'Gene', 'NormDepth']] # Kepp only columns of interest
data = data.loc[data['Gene'].isin(location['Gene']),:].copy() # Remove genes with no location information (mito genes)
data = data.loc[data['NormDepth'] > 0,].copy() # Remove genes with 0 Coverage

# Convert Genes and strain to integer IDs in order to optimize memory efficiency
strainsID = pd.DataFrame({'Strain': sorted(set(data['Strain']))})
strainsID['StrainID'] = strainsID.index + 1
data = data.merge(strainsID, on = 'Strain', how = 'left').drop(columns = ['Strain'])
genesID = pd.DataFrame({'Gene': sorted(set(data['Gene']))})
genesID['GeneID'] = genesID.index + 1
data = data.merge(genesID, on = 'Gene', how = 'left').drop(columns = ['Gene'])

# Process matrix
# ==============

# Decompose all possible variants
data['UpperCoverageBinEnd'] = data['NormDepth'].apply(lambda x: [y/2 + 0.25 for y in range(0, 1+2*(math.ceil(x)))])
#data['UpperCoverageBinEnd'] = data['NormDepth'].apply(lambda x: [y + 0.5 for y in range(0, 1+math.ceil(x))])
data = data.explode('UpperCoverageBinEnd')
data['UpperCoverageBinEnd'] = data['UpperCoverageBinEnd'].astype(float)

# Infer genotype
data['Genotype'] = 2
data.loc[data['NormDepth'] >= data['UpperCoverageBinEnd'], 'Genotype'] = 1

# Generate variant ID with haplotype info
variantsID = data.loc[:,['GeneID', 'UpperCoverageBinEnd']].drop_duplicates().reset_index(drop = True)
variantsID['VariantID'] = variantsID.index + 1
data = data.merge(variantsID, on = ['GeneID', 'UpperCoverageBinEnd'], how = 'left').drop(columns = ['GeneID', 'UpperCoverageBinEnd'])

# Add location information and order variants
variantsID = variantsID.merge(genesID, on = 'GeneID', how = 'left').merge(location, on = 'Gene', how = 'left')
variantsID = variantsID.sort_values(by = ['Chromosome', 'Start', 'Gene', 'UpperCoverageBinEnd'])

# Convert to ped
# ==============
# Pivot wider and add null (2) genotypes
ped = data.pivot(index = 'StrainID', columns = 'VariantID', values = 'Genotype').fillna(2).astype(int)
# Order variants and duplicate to simulate 2 haplotypes
ped = ped.loc[:,[x for x in variantsID['VariantID'] for _ in (0, 1)]]
# Add header columns
ped.index = [str(strainsID.loc[strainsID['StrainID'] == x, 'Strain'].iloc[0]) for x in ped.index]
ped.insert(loc = 0, column = 'Family', value = [x.split('_')[0] for x in ped.index])
ped.insert(loc = 1, column = 'Sample', value = [x.split('_')[1] if '_' in x else x for x in ped.index])
ped.insert(loc = 2, column = 'Paternal', value = 0)
ped.insert(loc = 3, column = 'Maternal', value = 0)
ped.insert(loc = 4, column = 'Sex', value = 0)
ped.insert(loc = 5, column = 'Affection', value = -9)

# Generate map file
# =================

map = variantsID.loc[:,['Gene', 'UpperCoverageBinEnd', 'Chromosome', 'Start']]
map['ID'] = map['Gene'] + '_Cov' + map['UpperCoverageBinEnd'].astype(str)
map['GeneticDistance'] = 0
map = map.loc[:,['Chromosome', 'ID', 'GeneticDistance', 'Start']]

# Write output files
# ==================
ped.to_csv(f'{outputPath}.ped', sep = '\t', index = False, header = False, quoting = csv.QUOTE_NONE)
map.to_csv(f'{outputPath}.map', sep = '\t', index = False, header = False, quoting = csv.QUOTE_NONE)

