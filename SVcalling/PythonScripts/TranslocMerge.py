#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/03/06
# version ='1.1'
# ---------------------------------------------------------------------------
'''
This script merges translocations together and output a Jasmine-like VCF.
'''
# ---------------------------------------------------------------------------
import pandas as pd
import argparse
import os
import itertools
import networkx as nx
from datetime import datetime
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
parser.add_argument("-v", "--VCF", help="Path of the single-sample VCFs", required=True, nargs = '+')
parser.add_argument("-kb", "--kb_margin", help="Maximum length between 2 breakpoints to merge translocations (kb)", required=True, type = int, default = 10)
parser.add_argument("-o", "--output", help="Output file (VCF format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

VCFPaths=[os.path.abspath(vcf) for vcf in args.VCF]
dist_kb=args.kb_margin
dist_bp = dist_kb * 1000
outputPath=os.path.abspath(args.output)

#VCFPaths = ['01_SV_VCF/AAA_output/AAA.SVs_all.IDset.vcf', '01_SV_VCF/AAB_2_output/AAB_2.SVs_all.IDset.vcf', '01_SV_VCF/AAB_output/AAB.SVs_all.IDset.vcf', '01_SV_VCF/AAC.HP1_output/AAC.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAC.HP2_output/AAC.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAD_2_output/AAD_2.SVs_all.IDset.vcf', '01_SV_VCF/AAD_output/AAD.SVs_all.IDset.vcf', '01_SV_VCF/AAE_2_output/AAE_2.SVs_all.IDset.vcf', '01_SV_VCF/AAE_output/AAE.SVs_all.IDset.vcf', '01_SV_VCF/AAI_2_output/AAI_2.SVs_all.IDset.vcf', '01_SV_VCF/AAI_output/AAI.SVs_all.IDset.vcf', '01_SV_VCF/AAL_2.HP1_output/AAL_2.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAL_2.HP2_output/AAL_2.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAL.HP1_output/AAL.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAL.HP2_output/AAL.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAM_2_output/AAM_2.SVs_all.IDset.vcf', '01_SV_VCF/AAM_output/AAM.SVs_all.IDset.vcf', '01_SV_VCF/AAN_2.HP1_output/AAN_2.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAN_2.HP2_output/AAN_2.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAN.HP1_output/AAN.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAN.HP2_output/AAN.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAP_2.HP1_output/AAP_2.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAP_2.HP2_output/AAP_2.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAP.HP1_output/AAP.HP1.SVs_all.IDset.vcf', '01_SV_VCF/AAP.HP2_output/AAP.HP2.SVs_all.IDset.vcf', '01_SV_VCF/AAQ_2_output/AAQ_2.SVs_all.IDset.vcf', '01_SV_VCF/AAQ_output/AAQ.SVs_all.IDset.vcf', '01_SV_VCF/AAR_2_output/AAR_2.SVs_all.IDset.vcf', '01_SV_VCF/AAR_output/AAR.SVs_all.IDset.vcf', '01_SV_VCF/AAS_output/AAS.SVs_all.IDset.vcf', '01_SV_VCF/AAT_2_output/AAT_2.SVs_all.IDset.vcf', '01_SV_VCF/AAT_output/AAT.SVs_all.IDset.vcf', '01_SV_VCF/AAV_2_output/AAV_2.SVs_all.IDset.vcf', '01_SV_VCF/AAV_output/AAV.SVs_all.IDset.vcf']

# Read input files ==========================================================
logging.info('Gathering translocations from all VCFs')
# For each VCF, only translocations are kept
translocations = pd.DataFrame()
totalStrains = []
for path in VCFPaths:
	strain = path.split('/')[-1].replace('.SVs_all.IDset.vcf', '')
	totalStrains += [strain]
	VCF = pd.read_csv(path, sep = '\t', comment = "#", names = ["CHROM1", "POS1", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
	# Filter translocations and get breakpoints
	VCF = VCF.loc[ VCF['ID'].str.contains('_TRA_'),:]
	if VCF.shape[0] > 0:
		VCF[['CHROM2', 'POS2']] = VCF['ALT'].str.replace('[', '').str.replace(']', '').str.split(':', expand = True)
		# Order breaks points
		VCF['SmallerBreakPoint'] = VCF.loc[:,['CHROM1', 'CHROM2']].apply(lambda x: x.str.replace('chromosome', '').astype(int).idxmin(), axis = 1)
		VCF.loc[VCF['CHROM1'] == VCF['CHROM2'],'SmallerBreakPoint'] = VCF.loc[:,['POS1', 'POS2']].apply(lambda x: x.astype(int).idxmin().replace('POS', 'CHROM'), axis = 1)
		## BreakPoint 1
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'BP1_CHR'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'CHROM1']
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'BP1_POS'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'POS1']
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'BP1_CHR'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'CHROM2']
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'BP1_POS'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'POS2']
		## BreakPoint 2
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'BP2_CHR'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'CHROM2']
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'BP2_POS'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM1', 'POS2']
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'BP2_CHR'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'CHROM1']
		VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'BP2_POS'] = VCF.loc[VCF['SmallerBreakPoint'] == 'CHROM2', 'POS1']
		# Concat
		translocations = pd.concat([translocations, VCF.loc[:,['ID', 'BP1_CHR', 'BP1_POS', 'BP2_CHR', 'BP2_POS']]])

translocations.index = translocations['ID']
translocations['BP1_POS'] = translocations['BP1_POS'].astype(int)
translocations['BP2_POS'] = translocations['BP2_POS'].astype(int)

# Add length
# The length of the translocation is defines as the distance between a breakpoint and 
# the closest extremity of the chromosome. Since there are 2 breakpoints, the larger
# length is kept. 
refChr = dict({'chromosome1' : 230218, 'chromosome2' : 813184, 'chromosome3' : 316620, 'chromosome4' : 1531933, 'chromosome5' : 576874, 'chromosome6' : 270161, 'chromosome7' : 1090940, 'chromosome8' : 562643, 'chromosome9' : 439888, 'chromosome10' : 745751, 'chromosome11' : 666816, 'chromosome12' : 1078177, 'chromosome13' : 924431, 'chromosome14' : 784333, 'chromosome15' : 1091291, 'chromosome16' : 948066})

translocations.loc[:,'BP1_CHR_LEN'] = translocations['BP1_CHR'].apply(lambda x: refChr.get(x))
translocations.loc[:,'BP1_DIST_END'] = translocations['BP1_CHR_LEN'] - translocations['BP1_POS']
translocations.loc[:,'BP1_LEN'] = translocations.loc[:,['BP1_POS', 'BP1_DIST_END']].min(axis = 1)

translocations.loc[:,'BP2_CHR_LEN'] = translocations['BP2_CHR'].apply(lambda x: refChr.get(x))
translocations.loc[:,'BP2_DIST_END'] = translocations['BP2_CHR_LEN'] - translocations['BP2_POS']
translocations.loc[:,'BP2_LEN'] = translocations.loc[:,['BP2_POS', 'BP2_DIST_END']].min(axis = 1)

translocations.loc[:,'LEN'] = translocations.loc[:,['BP1_LEN', 'BP2_LEN']].max(axis = 1)

# Translocation network =====================================================
# Translocations are linked if the breakpoints are on the same chromosome and
# positions are less than X kb away

logging.info("Build network")

# Compare each combination of translocations
combinations = pd.DataFrame(list(itertools.combinations(translocations['ID'], 2)), columns = ['ID1', 'ID2'])
combinations[['1_BP1_CHR', '1_BP1_POS', '1_BP2_CHR', '1_BP2_POS']] = translocations.loc[combinations['ID1'],['BP1_CHR', 'BP1_POS', 'BP2_CHR', 'BP2_POS']].reset_index(drop = True)
combinations[['2_BP1_CHR', '2_BP1_POS', '2_BP2_CHR', '2_BP2_POS']] = translocations.loc[combinations['ID2'],['BP1_CHR', 'BP1_POS', 'BP2_CHR', 'BP2_POS']].reset_index(drop = True)
# Remove combinations with different chromosomes
combinations = combinations.loc[(combinations['1_BP1_CHR'] == combinations['2_BP1_CHR']) & (combinations['1_BP2_CHR'] == combinations['2_BP2_CHR']),:]
#combinations.to_csv(outputPath+'_Combinations.csv')
# Remove combinations with breakpoints too distant
combinations = combinations.loc[(abs(combinations['1_BP1_POS'] - combinations['2_BP1_POS']) < dist_bp) & (abs(combinations['1_BP2_POS'] - combinations['2_BP2_POS']) < dist_bp),:]

# Create network
G = nx.Graph()
G.add_nodes_from(translocations['ID'].tolist())
G.add_edges_from([(combinations.loc[index,'ID1'], combinations.loc[index,'ID2']) for index in combinations.index])

# Format output VCF =========================================================
# Since a translocation does not have starting and ending point, the start is
# the first breakpoint and the end is the second breakpoint. 

logging.info("Create final VCF")

finalVCF = pd.DataFrame(columns = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
for c in list(nx.connected_components(G)):
	IDs = list(c)
	strains = [x.split(':')[0] for x in IDs]
	ID = IDs[0]
	# SV statistics
	CHR1=translocations.loc[IDs[0],'BP1_CHR']
	CHR2=translocations.loc[IDs[0],'BP2_CHR']
	AVG_START = int(translocations.loc[IDs,'BP1_POS'].mean().round())
	STARTVARIANCE = translocations.loc[IDs,'BP1_POS'].std()
	AVG_END = int(translocations.loc[IDs,'BP2_POS'].mean().round())
	ENDVARIANCE = translocations.loc[IDs,'BP2_POS'].std()
	AVG_LEN = int(translocations.loc[IDs,'LEN'].mean().round())
	# SV genotypes
	genotypes = pd.Series(0, index = totalStrains)
	genotypes.loc[strains] = 1
	genotypes = ''.join(genotypes.astype(str).tolist())
	# Build INFO field
	INFO = f'SVLEN={AVG_LEN};SVTYPE=TRA;STARTVARIANCE={STARTVARIANCE};ENDVARIANCE={ENDVARIANCE};AVG_LEN={AVG_LEN};AVG_START={AVG_START};AVG_END={AVG_END};SUPP_VEC_EXT={genotypes};IDLIST_EXT={",".join(IDs)};SUPP_EXT={len(IDs)};SUPP_VEC={genotypes};SUPP={len(IDs)};SVMETHOD=JASMINE_Homemade;IDLIST={",".join(IDs)};BREAKPOINT1={CHR1}:{AVG_START};BREAKPOINT2={CHR2}:{AVG_END}'
	newRecord = pd.DataFrame({"#CHROM": [CHR1], 
								"POS": [AVG_START], 
								"ID": [ID],
								"REF": ["."], 
								"ALT": ["<TRA>"], 
								"QUAL": ["."], 
								"FILTER": ["PASS"], 
								"INFO": [INFO], 
								"FORMAT": ["GT"], 
								"SAMPLE": ["1/1"]})
	finalVCF = pd.concat([finalVCF, newRecord])

# Order VCF
finalVCF['CNUM'] = finalVCF['#CHROM'].str.replace('chromosome', '').astype(int)
finalVCF = finalVCF.sort_values(['CNUM', 'POS'], ascending=[True, True]).reset_index(drop = True)
finalVCF = finalVCF.drop(columns = ['CNUM'])

# Add index number to SV IDs
finalVCF = finalVCF.reset_index()
finalVCF['ID'] = finalVCF['index'].astype(str) + '_' + finalVCF['ID']
finalVCF = finalVCF.drop(columns = ['index'])

# Get VCF header
header = []
sample = VCFPaths[0].split('/')[-1].replace('.SVs_all.IDset.vcf', '')
with open(VCFPaths[0], 'r') as file:
	for line in file:
		if line.startswith('##'):
			header += [line]
# Add Info fields
lastInfoIndex = [i for i, e in enumerate(header) if e.startswith("##INFO")][-1]
topHeader = header[:lastInfoIndex + 1]
bottomHeader = header[lastInfoIndex + 1:]
topHeader += ['##INFO=<ID=SUPP_VEC,Number=1,Type=String,Description="Vector of supporting samples">\n']
topHeader += ['##INFO=<ID=SUPP_VEC_EXT,Number=1,Type=String,Description="Vector of supporting samples, potentially extended across multiple merges">\n']
topHeader += ['##INFO=<ID=SUPP,Number=1,Type=Integer,Description="Number of samples supporting the variant">\n']
topHeader += ['##INFO=<ID=SUPP_EXT,Number=1,Type=Integer,Description="Number of samples supporting the variant, potentially extended across multiple merges">\n']
topHeader += ['##INFO=<ID=IDLIST,Number=.,Type=String,Description="Variant IDs of variants merged to make this call (at most 1 per sample)">\n']
topHeader += ['##INFO=<ID=IDLIST_EXT,Number=.,Type=String,Description="Variant IDs of variants merged, potentially extended across multiple merges">\n']
topHeader += ['##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="">\n']
topHeader += ['##INFO=<ID=STARTVARIANCE,Number=1,Type=String,Description="Variance of start position for variants merged into this one">\n']
topHeader += ['##INFO=<ID=ENDVARIANCE,Number=1,Type=String,Description="Variance of end position for variants merged into this one">\n']
topHeader += ['##INFO=<ID=AVG_START,Number=1,Type=String,Description="Average start position for variants merged into this one">\n']
topHeader += ['##INFO=<ID=AVG_END,Number=1,Type=String,Description="Average end position for variants merged into this one">\n']
topHeader += ['##INFO=<ID=AVG_LEN,Number=1,Type=String,Description="Average length for variants merged into this one">\n']
topHeader += ['##INFO=<ID=PRECISE,Number=0,Type=Flag,Description="Precise structural variation">\n']
topHeader += ['##INFO=<ID=IMPRECISE,Number=0,Type=Flag,Description="Imprecise structural variation">\n']
topHeader += ['##INFO=<ID=BREAKPOINT1,Number=1,Type=String,Description="First breakpoint of a translocation event">\n']
topHeader += ['##INFO=<ID=BREAKPOINT2,Number=1,Type=String,Description="Second breakpoint of a translocation event">\n']
header = topHeader + bottomHeader
# Add additionnal commands
header += [f'##CustomCommands=Translocation merging from {len(VCFPaths)} VCF files with {dist_kb}kb breakpoint margin; Date = {datetime.now().strftime("%a %b %d %H:%M:%S %Y")}\n']

# Write final file ==========================================================
with open(outputPath, 'w') as out:
	for line in header:
		out.write(line)

finalVCF = finalVCF.rename(columns = {'SAMPLE': sample})
finalVCF.to_csv(outputPath, sep = '\t', mode='a', index = False)




