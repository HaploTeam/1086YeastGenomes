#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/03/06
# version ='1.2'
# ---------------------------------------------------------------------------
'''
This script transfers annotation from a reference to a draft assembly using 
blastn. For each gene of the draft assembly, the script looks for an ortholog
in the reference genes. Orthology is indicated if more than 80% of the 
reference gene is present with >= 95% identity in the draft gene. If a smaller
alignement is found, the gene is annotated as truncated. 

It takes as argument: 
- the reference CDS
- the prefix of gff3 and fasta (cds, trimmed_cds and pep) from LRSDAY
- the output prefix of gff3 and fasta
'''
# ---------------------------------------------------------------------------
import os
import subprocess
import argparse
import random
import pandas as pd
import sys
sys.path.insert(1, '/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
# ---------------------------------------------------------------------------
# Definitions
def addInfo(dataFrame):
	"""This function uses a blast alignment input and adds information about
	coverage (%) of the query sequence by the subject sequence and opposite, 
	as well as information of alignment spanning the beginning and the end of 
	the subject sequence. """
	# Compute query coverage ================================================
	BEDgeneQ = BED([BEDcoordinates('Gene', 0, dataFrame.iloc[0,dataFrame.columns.get_loc('qlen')])])
	BEDalignmentsQ = []
	for i in dataFrame.index:
		BEDalignmentsQ += [BEDcoordinates('Gene', dataFrame.loc[i,'qstart']-1, dataFrame.loc[i,'qend'])]
	BEDalignmentsQ = BED(BEDalignmentsQ)
	qcov = BEDgeneQ.overlapLen(BEDalignmentsQ, percent = True)
	# Compute subject coverage ==============================================
	BEDgeneS = BED([BEDcoordinates('Gene', 0, dataFrame.iloc[0,dataFrame.columns.get_loc('slen')])])
	BEDalignmentsS = []
	for i in dataFrame.index:
		BEDalignmentsS += [BEDcoordinates('Gene', dataFrame.loc[i,'sstart']-1, dataFrame.loc[i,'send'])]
	BEDalignmentsS = BED(BEDalignmentsS)
	scov = BEDgeneS.overlapLen(BEDalignmentsS, percent = True)
	# Check if alignments span begining of subject ==========================
	firstPos = min(pd.concat([dataFrame.loc[:,'sstart'], dataFrame.loc[:,'send']]))
	lastPos = max(pd.concat([dataFrame.loc[:,'sstart'], dataFrame.loc[:,'send']]))
	startSubject = True if firstPos <= 10 else False
	endSubject = True if lastPos >= dataFrame.iloc[0,dataFrame.columns.get_loc('slen')]-10 else False
	# Retrieve minimu evalue ================================================
	evalue = min(dataFrame.loc[:,'evalue'])
	# Return single line df with all informations
	return pd.DataFrame({'evalue': [evalue], 'scov': [scov], 'qcov': [qcov], 'startSubject': [startSubject], 'endSubject': [endSubject]})

def findGeneCorrespondance(dataFrame):
	# If some subject are complete, juste keep complete genes
	if 'complete' in dataFrame['suffix'].unique():
		dataFrame = dataFrame.loc[dataFrame['suffix'] == 'complete',:]
		# Keep genes with minimum evalue and maximum query coverage
		dataFrame = dataFrame.loc[dataFrame['evalue'] == min(dataFrame.loc[:,'evalue']),:]
		dataFrame = dataFrame.loc[dataFrame['qcov'] == max(dataFrame.loc[:,'qcov']),:]
		# Final gene corespondance
		dataFrame = dataFrame.sort_values(by = 'sseqid')
		annot = "/".join(dataFrame['sseqid'].tolist())
	# If no complete gene are present, but truncated genes are
	elif 'truncated' in dataFrame['suffix'].unique():
		dataFrame = dataFrame.loc[dataFrame['suffix'] == 'truncated',:]
		# Keep genes with minimum evalue and maximum query coverage
		dataFrame = dataFrame.loc[dataFrame['evalue'] == min(dataFrame.loc[:,'evalue']),:]
		dataFrame = dataFrame.loc[dataFrame['qcov'] == max(dataFrame.loc[:,'qcov']),:]
		# Final gene corespondance
		dataFrame = dataFrame.sort_values(by = 'sseqid')
		dataFrame['sseqid'] = dataFrame['sseqid'] + '_truncated'
		annot = "/".join(dataFrame['sseqid'].tolist())
	# If no complete nor truncated genes are present
	else:
		annot = dataFrame.iloc[0,dataFrame.columns.get_loc('qseqid')]
	return pd.DataFrame({'annot': [annot]})

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-r", "--referenceORFs", help="Reference ORFs fasta file", required=True)
parser.add_argument("-f", "--files", help="Files prefix from LRSDAY output (fasta and gff3)", required=True)
parser.add_argument("-o", "--output", help="Output prefix for fasta and gff3", required=True)

# Read arguments from the command line
args = parser.parse_args()

referenceORFs=os.path.abspath(args.referenceORFs)
files=os.path.abspath(args.files)
output=os.path.abspath(args.output)

#referenceORFs='/home/vloegler/Pangenome3/03-data/orf_genomic_all_R64-3-1_20210421.fasta'
#files='/home/vloegler/Pangenome3/03-data/CollapsedGenomesAnnotations/AAD.nuclear_genome.Final'
#output='/home/vloegler/Pangenome3/04-analysis/test'

# Create temporary workdir
prefix = files.split("/")[-1]
workdir = f'Workdir_TransferAnnotations_{prefix}_{random.randint(10000, 99999)}'
os.mkdir(workdir)
os.chdir(workdir)

# Make blast db with reference ORFs
referenceLink = 'ReferenceORFs.fasta'
os.symlink(referenceORFs, referenceLink)
command = f'makeblastdb -dbtype nucl -in {referenceLink} -title {referenceLink}'.split()
subprocess.call(command)

# Copy cds from draft and keep only first ID for each gene
CDS = Fasta(f'{files}.cds.fa')
for seq in CDS:
	seq.id = seq.id.split("|")[0]
	seq.description = f'>{seq.id}\n'
CDS.toFile('Query.fasta')

# Run blastn between draft trimmed_CDS and reference ORFs
command = ['blastn', 
	   '-db', 'ReferenceORFs.fasta', 
	   '-query', 'Query.fasta', 
	   '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen', 
	   '-num_threads', '8', 
	   '-dust', 'no', 
	   '-perc_identity', '95', 
	   '-strand', 'plus', 
	   '-out', 'Output.blastn']
subprocess.call(command)

# Read blast output
blastOut = pd.read_csv('Output.blastn', sep = '\t', header = None, names = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'.split())

# Group by query subject pairs and add coverage and location information
blastOutGrouped = blastOut.groupby(['qseqid', 'sseqid', 'qlen', 'slen']).apply(addInfo).reset_index()

# For each query subject pair, check if subject sequence is complete, truncated or not present
# Annotate as truncated when subject coverage is 50%
blastOutGrouped.loc[(blastOutGrouped['scov'] >= 50),'suffix'] = 'truncated'
# Annotate as truncated when total query coverage is larger than 90%
blastOutGrouped.loc[(blastOutGrouped['qcov'] >= 90),'suffix'] = 'truncated'
# Annotate as complete when beginning AND end of subject is covered, and total coverage is 30%
blastOutGrouped.loc[(blastOutGrouped['scov'] >= 30) & blastOutGrouped['startSubject'] & blastOutGrouped['endSubject'],'suffix'] = 'complete'
# Annotate as complete when beginning AND end of subject is covered, and more than 80% of query is covered
blastOutGrouped.loc[(blastOutGrouped['qcov'] >= 80) & blastOutGrouped['startSubject'] & blastOutGrouped['endSubject'],'suffix'] = 'complete'
# Annotate as complete when subject sequenced is covered at 80%
blastOutGrouped.loc[blastOutGrouped['scov'] >= 80,'suffix'] = 'complete'

# Find reference gene correspondance
blastOutGrouped = blastOutGrouped.groupby(['qseqid']).apply(findGeneCorrespondance).reset_index()

# Add genes with no blast results
missingGenes = [x for x in CDS.getID() if not x in blastOutGrouped.loc[:,'qseqid'].tolist()]
finalAnnotation = pd.concat([blastOutGrouped, pd.DataFrame({'qseqid': missingGenes, 'annot': missingGenes})]).reset_index(drop = True)

# convert to dictionnary
finalAnnotation = finalAnnotation.set_index(finalAnnotation['qseqid']).loc[:,'annot'].to_dict()

# Change names in files
## cds fasta file
cds = Fasta(f'{files}.cds.fa')
for seq in cds:
	IDlist = seq.id.split("|")
	IDlist[2] = finalAnnotation[IDlist[0]]
	seq.id = '|'.join(IDlist)
	seq.description = f'>{seq.id}\n'
cds.toFile(f'{output}.cds.fa')
print(f'CDS fasta written to {output}.cds.fa')

## trimmed cds fasta file
try:
	trimmedcds = Fasta(f'{files}.trimmed_cds.fa')
	for seq in trimmedcds:
		IDlist = seq.id.split("|")
		IDlist[2] = finalAnnotation[IDlist[0]]
		seq.id = '|'.join(IDlist)
		seq.description = f'>{seq.id}\n'
	trimmedcds.toFile(f'{output}.trimmed_cds.fa')
	print(f'Trimmed CDS fasta written to {output}.trimmed_cds.fa')
except:
	pass

## pep fasta file
pep = Fasta(f'{files}.pep.fa')
for seq in pep:
	IDlist = seq.id.split("|")
	IDlist[2] = finalAnnotation[IDlist[0]]
	seq.id = '|'.join(IDlist)
	seq.description = f'>{seq.id}\n'
pep.toFile(f'{output}.pep.fa')
print(f'Pep fasta written to {output}.pep.fa')

## gff3 file
with open(f'{files}.gff3', 'r') as gff3:
	with open(f'{output}.gff3', 'w') as out:
		for line in gff3:
			if len(line.split()) > 2 and line.split()[2] == 'gene':
				line = line.split()
				ID = line[8].split(";")[0].split("=")[1]
				line[8] = f'ID={ID};Name={finalAnnotation[ID]}'
				line = '\t'.join(line)+'\n'
			out.write(line)
print(f'Annotation GFF3 file written to {output}.gff3')

