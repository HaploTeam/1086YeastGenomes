#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2024/01/23
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script takes as input GWAS results files and genotypes and compute 
linkage groups. Linkage groups are computed across all significant variants, 
and can be additionnally clustered by type of variants with option -p. 

Linkage groups are formed using a graph-based approach using plink LD 
results. 
'''
# ---------------------------------------------------------------------------
import re
import os
import subprocess
import pandas as pd
import networkx as nx
import argparse
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = '''
This script takes as input GWAS results files and genotypes and compute 
linkage groups. Linkage groups are computed across all significant variants, 
and can be additionnally clustered by type of variants with option -p. 

Linkage groups are formed using a graph-based approach using plink LD 
results. 
''')
parser.add_argument("-g", "--genotype", help="Genotype matrix in Plink format", required=True)
parser.add_argument("-r", "--results", help="GWAS results (Pheno.signif_snps.txt file)", required=True)
parser.add_argument("-w", "--window_size", help="Maximum distance between 2 variants to be considered in a same LD group (in kb)", type = int, default = 50)
parser.add_argument("-r2", "--r2", help="LD threshold above which 2 variants will be put in the same LD group", type = float, default = 0.5)
parser.add_argument("-p", "--prefix_type", help="Prefix of different type of variants. If precised (ex: --prefix_type SNP INDEL SV CNV), linkage groups will be built independentely on each type of variant, by discriminating them over the prefix of the variant ID. In addition, linkage groups acroos all significant variant will still be added. ", required=False, nargs='*')

# Read arguments from the command line
args = parser.parse_args()
genotypePath = os.path.abspath(args.genotype)
resultsPath = os.path.abspath(args.results)
windowSize = args.window_size
r2threshold = args.r2
variantTypes = args.prefix_type

# ====================================
# Reads GWAS results
results = pd.read_csv(resultsPath, sep = '\t', header = 0)

if results.shape[0] > 1:
	# Extract significant variants and get linkage
	results.to_csv(resultsPath + '.OnlyVariantID.txt', header = False, index = False, columns = ['SNP'])
	cmd = ['plink', 
		'--bfile', genotypePath, 
		'--extract', resultsPath+'.OnlyVariantID.txt', 
		'--out', resultsPath+'.OnlyVariantID.plink', 
		'--r2', 
		'--ld-window-r2', str(r2threshold), 
		'--ld-window', str(1000000), 
		'--ld-window-kb', str(windowSize)]
	subprocess.call(cmd)
	linkage = pd.read_csv(resultsPath+'.OnlyVariantID.plink.ld', header = 0, delim_whitespace=True)
	cmd = ['rm', '-f', resultsPath+'.OnlyVariantID.txt', resultsPath+'.OnlyVariantID.plink.ld', resultsPath+'.OnlyVariantID.plink.log', resultsPath+'.OnlyVariantID.plink.nosex']
	subprocess.call(cmd)

	# Create linkage groups
	G = nx.Graph()
	G.add_nodes_from(results['SNP'])
	edges = [(linkage.loc[i,'SNP_A'], linkage.loc[i,'SNP_B']) for i in linkage.index]
	G.add_edges_from(edges)

	# Add cluster information to the results table
	clusters = list(nx.connected_components(G))
	results['LinkageGroup'] = 0
	i = 1
	for c in clusters:
		results.loc[results['SNP'].isin(c), 'LinkageGroup'] = f'LD_{i}'
		i += 1

	# If variant prefixes are precised, compute linkage groups for each type of variants
	if not variantTypes is None:
		results['LinkageGroup_PerVariantType'] = 0
		for type in variantTypes:
			# Subset linkage per type
			linkage_sub = linkage.loc[linkage['SNP_A'].str.startswith(type) & linkage['SNP_B'].str.startswith(type),:]
			# Create linkage groups
			G = nx.Graph()
			G.add_nodes_from(results.loc[results['SNP'].str.startswith(type),'SNP'])
			edges = [(linkage_sub.loc[i,'SNP_A'], linkage_sub.loc[i,'SNP_B']) for i in linkage_sub.index]
			G.add_edges_from(edges)
			# Add cluster information to the results table
			clusters = list(nx.connected_components(G))
			i = 1
			for c in clusters:
				results.loc[results['SNP'].isin(c), 'LinkageGroup_PerVariantType'] = f'LD_{type}_{i}'
				i += 1
elif results.shape[0] == 1:
	results['LinkageGroup'] = 'LD_1'
	if len(variantTypes) > 0:
		for type in variantTypes:
			results.loc[results['SNP'].str.startswith(type), 'LinkageGroup_PerVariantType'] = f'LD_{type}_1'
else:
	results['LinkageGroup'] = None
	if len(variantTypes) > 0:
		results['LinkageGroup_PerVariantType'] = None

results.to_csv(f'{re.sub(".txt", "", resultsPath)}.LinkageGroups.txt', sep = '\t', index = False)
