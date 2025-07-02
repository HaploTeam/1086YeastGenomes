#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/03/06
# version ='1.1'
# ---------------------------------------------------------------------------
'''
This script converts the output of Jasmine (merging of single-sample SV VCFs)
into a standardized multi-sample VCF. 
It also retrieves translocation breakpoints. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import re
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
parser.add_argument("-v", "--jasmineVCF", help="Path of the VCF obtained with Jasmine", required=True)
parser.add_argument("-i", "--jasmineInput", help="Input file given to Jasmine", required=True)
parser.add_argument("-d", "--singleVCFdir", help="Directory containing all the single sample VCF used for the Jasmine merging", required=True)
parser.add_argument("-o", "--output", help="Output file (VCF format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

jasmineVCFPath=os.path.abspath(args.jasmineVCF)
jasmineInputPath=os.path.abspath(args.jasmineInput)
singleVCFdir=os.path.abspath(args.singleVCFdir)
outputPath=os.path.abspath(args.output)

#jasmineVCFPath=os.path.abspath("../02_JasmineMerging/SV.JasmineMerged.vcf")
#jasmineInputPath=os.path.abspath("../02_JasmineMerging/listVCFs.txt")

# Read input files ==========================================================
# Read VCF # Coluimn ALT is replaced by ALT_ALLELE because ALT is a strain ID
VCF = pd.read_csv(jasmineVCFPath, sep = '\t', comment = "#", names = ["#CHROM", "POS", "ID", "REF", "ALT_ALLELE", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
VCF = VCF.loc[:,["#CHROM", "POS", "ID", "REF", "ALT_ALLELE", "QUAL", "FILTER", "INFO", "FORMAT"]] # Remove SAMPLE column
# Reads strains
strains = pd.read_csv(jasmineInputPath, header = None).iloc[:,0].tolist()
strains = [re.sub('.SVs_all.IDset.vcf', '', x) for x in strains]
strainsBest = [x for x in strains if not x.split('.')[0].endswith('_2')]

# Extract genotype for each strain ==========================================
def getGenotypes(INFO):
	SUPP_VEC = [x for x in INFO.split(';') if x.startswith('SUPP_VEC=')][0].split("=")[1] # Get SUPP_VEC INFO field
	return pd.Series([int(GT) for GT in SUPP_VEC], index = strains)
VCF.loc[:,strains] = VCF['INFO'].apply(getGenotypes)
VCF = VCF.copy() # de-fragmentate dataframe

# Remove SV only present in second best assemblies ==========================
VCF.loc[:,'SUM_BEST'] = VCF.loc[:,strainsBest].apply(sum, axis = 1)
VCF = VCF.loc[VCF['SUM_BEST'] > 0,:]

# Remove singletons that are absent from the second best assembly ===========
VCF.loc[:,'FalsePos'] = 0
for s in strainsBest:
	# Second best strain
	s2 = '_2.'.join(s.split('.')) if '.HP' in s else f'{s}_2'
	if s2 in VCF.columns:
		# Add false positive information: variant is singleton, present in best assembly and absent from second best assembly
		VCF.loc[(VCF['SUM_BEST'] == 1) & (VCF[s] == 1) & (VCF[s2] == 0),'FalsePos'] = 1
	else: # If no Second Best assembly, all singletons for this strain are considered as false positives
		VCF.loc[(VCF['SUM_BEST'] == 1) & (VCF[s] == 1),'FalsePos'] = 1
# Remove False positives
VCF = VCF.loc[VCF['FalsePos'] == 0,:]

# Format final VCF ==========================================================
# Get phased genotypes
strainsUnique = sorted(set(sorted([re.sub('.HP2', '', re.sub('.HP1', '', x)) for x in strainsBest])))
finalVCF = pd.concat([VCF.loc[:,["#CHROM", "POS", "ID", "REF", "ALT_ALLELE", "QUAL", "FILTER", "INFO", "FORMAT"]], pd.DataFrame(index = VCF.index, columns = strainsUnique)], axis = 1)
for s in strainsUnique:
	if s in VCF.columns:
		finalVCF.loc[:,s] = VCF[s].astype(str) + '|' + VCF[s].astype(str)
	else:
		hp1 = f'{s}.HP1'
		hp2 = f'{s}.HP2'
		finalVCF.loc[:,s] = VCF[hp1].astype(str) + '|' + VCF[hp2].astype(str)

# Get VCF header
header = []
with open(jasmineVCFPath, 'r') as file:
	for line in file:
		if line.startswith('##'):
			header += [line]
# Add additionnal commands
vcfFiles = pd.read_csv(jasmineInputPath, header = None).iloc[:,0].tolist()
header += [f'##JasmineCommand=Merge {len(strains)} haplotypes from files {" ".join(vcfFiles)};\n']
header += [f'##CustomCommands=Convert Jasmine output to classic VCF file, using second best assemblies as confirmation for singletons;\n']

# Write final file ==========================================================
finalVCF = finalVCF.rename(columns={"ALT_ALLELE": "ALT"})
with open(outputPath, 'w') as out:
	for line in header:
		out.write(line)
finalVCF.to_csv(outputPath, sep = '\t', mode='a', index = False)
