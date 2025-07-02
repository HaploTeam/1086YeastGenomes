#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/06/28
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script select the 500 haplotypes containing the largest number of SVs. 
This 500 haplotypes will be used to build a graph pangenome using Minigraph 
Cactus. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import subprocess
import os
import shutil
import argparse
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = '''
This script select the 500 haplotypes containing the largest number of SVs. 
This 500 haplotypes will be used to build a graph pangenome using Minigraph 
Cactus. 
''')
parser.add_argument("-v", "--vcf", help="VCF file of the structural variants", required=True)
parser.add_argument("-s", "--strains", help="File with strains having high quality genome assemblies (one strain per line)", required=True)
parser.add_argument("-d", "--datadir", help="Path to the data directory", required=True)
parser.add_argument("-o", "--outdir", help="Path to the output directory. If existant, will be removed. ", required=True)


# Read arguments from the command line
args = parser.parse_args()

vcfPath = args.vcf
strainsPath = args.strains
datadir = args.datadir
outdir = args.outdir

# Import SV VCF file
header = subprocess.check_output(f'zgrep \#CHROM {vcfPath}'.split()).decode('utf-8').strip().split()
header[4] = 'ALT_ALLELE' # Change ALT column to ALT_ALLELE because there is an ALT strain and duplicate columns names are not allowed
strains = header[9:len(header)]
VCF = pd.read_csv(vcfPath, sep = "\t", comment = "#", names = header)
VCF = VCF.loc[:,strains] # Keep only genotype columns
TotalNbSVs = VCF.shape[0]

# Get high quality strains
HQ_strains = pd.read_csv(strainsPath, header = None).loc[:,0].tolist()

# Split haplotypes
for s in strains:
	if '1|0' in VCF[s].tolist() or '0|1' in VCF[s].tolist():
		newDF = pd.DataFrame(index = VCF.index)
		newDF[[f'{s}.HP1', f'{s}.HP2']] = VCF[s].str.split('|', expand=True)
		newDF = newDF.astype({f'{s}.HP1':"int",f'{s}.HP2':"int"})
		VCF = pd.concat([VCF.loc[:,[x for x in VCF.columns if x != s]], newDF], axis = 1)
	else:
		VCF.loc[VCF[s] == "0|0",s] = 0
		VCF.loc[VCF[s] == "1|1",s] = 1

# Remove DataFrame fragmentation
VCF = VCF.copy()

# Get haplotypes
haplotypes = VCF.columns.tolist()
HQ_haplotypes = [x for x in haplotypes if x.split('.')[0] in HQ_strains]

# Create directory for selected genomes
if os.path.exists(outdir):
	shutil.rmtree(outdir)

os.mkdir(outdir)

# Copy reference genome
subprocess.run(f'ln -s {datadir}/Sace_S288c_reference_FullMatrixID.fna {outdir}/'.split())

# First add high quality haplotypes
# While number of genome is below 500, take genomes with the largest number of new SVs

while len(os.listdir(outdir)) < 500 and VCF.shape[0] > 0:
	# Get the number of SV for each strain
	sumSV = VCF.loc[:,HQ_haplotypes].apply(sum, axis = 0).tolist()
	if max(sumSV) > 0: # If some high quality haplotype contains new SVs
		haploToAdd = HQ_haplotypes[sumSV.index(max(sumSV))] # Haplotype with the higher number of SVs
	else:
		sumSV = VCF.loc[:,haplotypes].apply(sum, axis = 0).tolist()
		haploToAdd = haplotypes[sumSV.index(max(sumSV))] # Haplotype with the higher number of SVs
	print(f'{haploToAdd}\t{max(sumSV)}')
	subprocess.run(f'ln -s {datadir}/{haploToAdd}.Final.fasta {outdir}/'.split())
	# Removes SV that are in this strain
	VCF = VCF.loc[VCF[haploToAdd] != 1,:]

print(f'{TotalNbSVs - VCF.shape[0]}/{TotalNbSVs} SVs added added with 500 haplotypes. ')
