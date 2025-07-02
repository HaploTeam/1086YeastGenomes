#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/03/06
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script filters the results of the blastp of non reference genes of Sace
pangenome on the RefSeq database. It infers the function of each protein 
based on sequence similarity. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import os
import csv
import sys
import argparse
import subprocess
import logging
logging.basicConfig(level=logging.INFO)
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
parser.add_argument("-f", "--pangenomeFasta", help="Fasta file of the pangenome to annotate", required=True)
parser.add_argument("-b", "--blastResults", help="Blast results", required=True)
parser.add_argument("-o", "--output", help="Output file (tsv format)", required=True)
parser.add_argument("-d", "--datasets", help="Path to the datasets exe file (https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/)", required=True)
parser.add_argument("-g", "--gene2go", help="Path to the gene2go.gz file containing correspondance between RefSeq GeneID and assigned GO terms (https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz)", required=True)

# Read arguments from the command line
args = parser.parse_args()

pangenomeFastaPath=os.path.abspath(args.pangenomeFasta) # /home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa
blastResultsPath=os.path.abspath(args.blastResults) # /home/vloegler/Pangenome_Sace1000ONT/04-analysis/10_FunctionInference/PangenomeNonRef.blastp
output=os.path.abspath(args.output)
datasets = os.path.abspath(args.datasets)
gene2goPath = os.path.abspath(args.gene2go)

#pangenomeFastaPath="/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa"
#blastResultsPath="/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/08_BlastOnRefSeq/PangenomeNonRef_on_RefSeq_Shen2018_Yue2017.blastp"
#datasets="/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts/datasets"
#gene2goPath="/shared/home/vloegler/Pangenome_Sace1000ONT/03-data/gene2go.gz"

# Get sequences IDs for the non ref genes of the pangenome
nonRefGenes = [x for x in Fasta(pangenomeFastaPath).getID() if x.startswith('YX')]

# Read blast results
blastResults = pd.read_csv(blastResultsPath, sep ='\t', header = None, names='qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames scomnames staxids'.split())
# Remove hits that are not in RefSeq
blastResults = blastResults.loc[blastResults['sseqid'].str.startswith('ref|'),:]

# Add subject coverage information
# ================================
# Convert Alignments and genes into BED objects (s for subject)
blastResults.loc[:,'sBED'] = blastResults['slen'].apply(lambda x: BED(BEDcoordinates('Gene', 0, x)))
blastResults.loc[:,'sA_BED'] = blastResults.loc[:,['sstart', 'send']].apply(lambda x: BED(BEDcoordinates('Gene', x.iloc[0]-1, x.iloc[1])), axis = 1)
# Group by QuerySubject pairs and sum alignment BEDs
blastResultsGrouped = blastResults.groupby(['qseqid', 'sseqid']).agg({'sBED': 'first', 'sA_BED':'sum'}).reset_index()
# Compute query and subject coverage
blastResultsGrouped.loc[:,'scovq'] = blastResultsGrouped.apply(lambda x: x.sBED.overlapLen(x.sA_BED, percent = True), axis = 1)
# Remove unecessary columns
blastResults = blastResults.drop(columns=['sBED', 'sA_BED'])
# Add subject coverage information
blastResults = pd.merge(blastResults, blastResultsGrouped.loc[:,['qseqid', 'sseqid', 'scovq']],  how='left', left_on=['qseqid','sseqid'], right_on = ['qseqid','sseqid'])

# Filter blast results
# ====================
minIdent=30
minCov=50

# Basic filters
blastResults = blastResults.loc[blastResults['pident'] >= minIdent,:]
#blastResults = blastResults.loc[blastResults['scovq'] >= minCov,:]
blastResults = blastResults.loc[blastResults['qcovs'] >= minCov,:]
blastResultsAll = blastResults # Keep all results for uncharaterized proteins information

# Remove hypothetical and uncharacterized proteins if not SGD records (Saccharomyces cerevisiae S288c)
blastResults = blastResults.loc[(blastResults['stitle'].str.lower().str.contains("hypothetical") == False) | (blastResults['sscinames'] == 'Saccharomyces cerevisiae S288C'),:]
blastResults = blastResults.loc[(blastResults['stitle'].str.lower().str.contains("uncharacterized") == False) | (blastResults['sscinames'] == 'Saccharomyces cerevisiae S288C'),:]

# Keep only first results for each protein
blastResults = blastResults.groupby(['qseqid']).first().reset_index()

# Try to get gene name from Uniprot
def getGene(ID):
	ID = ID.split('|')[1]
	logging.info(f'Searching Uniprot data for: {ID}')
	try: 
		uniprotOutput = subprocess.check_output(['curl', '-H', 'Accept: text/plain; format=flatfile', f'https://rest.uniprot.org/uniprotkb/search?query=(xref:refseq-{ID})']).decode(sys.stdout.encoding).split('\n')
	except subprocess.CalledProcessError:
		logging.info('Ran into subprocess.CalledProcessError, trying a second time')
		uniprotOutput = subprocess.check_output(['curl', '-H', 'Accept: text/plain; format=flatfile', f'https://rest.uniprot.org/uniprotkb/search?query=(xref:refseq-{ID})']).decode(sys.stdout.encoding).split('\n')
	if len(uniprotOutput) > 1:
		ID = [x for x in uniprotOutput if x.startswith('AC')][0].split()[1].split(';')[0]
		ProtName = [x for x in uniprotOutput if x.startswith('DE')]
		if len(ProtName) > 1:
			ProtName = [x for x in ProtName if 'RecName' in x or 'SubName' in x]
		ProtName = ProtName[0].split(';')[0].split('=')[1].split('{')[0].strip()
		GeneName = [x for x in uniprotOutput if x.startswith('GN')][0]
		GeneName1 = GeneName.split(';')[0].split('=')[1].split('{')[0].strip()
		if 'LocusName' in GeneName:
			GeneName2 = [x for x in GeneName.split(';') if 'OrderedLocusNames' in x][0].split('=')[1].split('{')[0].strip()
			if GeneName1 != GeneName2:
				GeneName1 = f'{GeneName1} {GeneName2}'
		return([ID, GeneName1, ProtName])
	else:
		return('')

blastResults['uniprotData'] = blastResults['sseqid'].apply(getGene)

# Final gene annotation =========================================================
# [Similar / Highly similar] to [database] [ID] [Species] [Gene name] [Prot name] [Function] [GO terms associated]

# [Similar / Highly similar]
blastResults['similarity'] = "Similar to"
blastResults.loc[blastResults['pident'] >= 80,'similarity'] = "Highly similar to"

# [database]
blastResults['db'] = "Uniprot"
blastResults.loc[blastResults['uniprotData'] == '', 'db'] = "RefSeq"

# [ID]
blastResults['ID'] = blastResults['sseqid'].str.replace('ref','').str.replace('|','')

def getFirst(x):
	return(x[0])
blastResults.loc[blastResults['db'] == 'Uniprot', 'ID'] = blastResults.loc[blastResults['db'] == 'Uniprot','uniprotData'].apply(getFirst)

# [Species]
def getSpecies(x):
	return(x.split(';')[0])
blastResults['Species'] = blastResults['sscinames'].apply(getSpecies)

# [names]
def getGene(ID): # Try to get gene data from RefSeq
	logging.info(f'Getting gene name for RefSeq protein: {ID}')
	ID = ID.split('|')[1]
	genbank = ''
	while genbank == '':
		esearch = subprocess.Popen(['esearch', '-db', 'protein', '-query', ID], stdout=subprocess.PIPE, text=True)
		efetch = subprocess.Popen(['efetch', '-format', 'genbank'], stdin=esearch.stdout, stdout=subprocess.PIPE, text=True)
		genbank, error = efetch.communicate()
	genbank = [x.strip() for x in genbank.split('\n')]
	gene = [x for x in genbank if x.startswith('/gene="')]
	locus_tag = [x for x in genbank if x.startswith('/locus_tag="')]
	toreturn = ''
	if len(gene) > 0:
		toreturn += gene[0].split("=")[1].strip('\"') + ' '
	if len(locus_tag) > 0:
		toreturn += locus_tag[0].split("=")[1].strip('\"') + ' '
	return(toreturn)
def getNames(x):
	return(x.split('[')[0].strip())
blastResults.loc[blastResults['db'] == 'RefSeq', 'names'] =  blastResults.loc[blastResults['db'] == 'RefSeq','sseqid'].apply(getGene) + blastResults.loc[blastResults['db'] == 'RefSeq','stitle'].apply(getNames)

def getNames(x):
	return(f'{x[1]} {x[2]}')
blastResults.loc[blastResults['db'] == 'Uniprot', 'names'] = blastResults.loc[blastResults['db'] == 'Uniprot','uniprotData'].apply(getNames)

# [Final]
blastResults['Final'] = blastResults['similarity'].astype(str) + ' ' + blastResults['db'].astype(str) + ' ' + blastResults['ID'].astype(str) + ' ' + blastResults['Species'].astype(str) + ' ' + blastResults['names'].astype(str)

# [GO terms associated]
gene2go = pd.read_csv(gene2goPath, sep = '\t') # Read file with gene_id / GO term equivalence
blastResults['RefSeqID'] = blastResults['sseqid'].str.replace('ref','').str.replace('|','')
def getGoTerms(AccessionNumber):
	ids = subprocess.check_output([datasets, "summary", "gene", "accession", AccessionNumber]).decode()
	if '"gene_id":' in ids:
		gene_id = int([x for x in ids.split(',') if '"gene_id":' in x][0].split(':')[1].replace('"', ''))
		logging.info(f'{AccessionNumber}\t{gene_id}')
		GOterms = gene2go.loc[gene2go['GeneID'] == gene_id, 'GO_ID']
		return(' '.join(set(GOterms)))
	else:
		logging.info('Skipped')
		return('')

blastResults['TransferredGOterms'] = blastResults['RefSeqID'].apply(getGoTerms)

# Find uncharacterized and dubious proteins
# Uncharacterized have a match in RefSeq but no function, while dubious have no match in RefSeq

blastResultsUnchara = blastResultsAll.loc[ ~blastResultsAll['qseqid'].isin(blastResults['qseqid'].tolist())]
UncharaGenes = set(sorted(blastResultsUnchara['qseqid'].tolist()))

noHit = pd.DataFrame(columns = blastResults.columns) 
noHit.loc[:,'qseqid'] = [x for x in nonRefGenes if not x in blastResults['qseqid'].tolist()]
noHit.loc[:,'Final'] = 'Dubious protein'
noHit.loc[:,'TransferredGOterms'] = ''
noHit.loc[noHit['qseqid'].isin(UncharaGenes),'Final'] = 'Uncharacterized protein'

functions = pd.concat([blastResults, noHit])
functions = functions.loc[:,['qseqid', 'Final', 'TransferredGOterms']]
functions.to_csv(output, sep = '\t', index = False, quoting = csv.QUOTE_NONNUMERIC)
