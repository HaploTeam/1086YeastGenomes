#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/07/26
# version ='1.1'
# ---------------------------------------------------------------------------
'''
This script remove some genes of the pangenome that show high level of 
exact redundancy with other genes of the pangenome. 
Exact redundancy is find using blastn. Gene are removed in an iterative 
manner. The script output processed nucleotide and protein fastas, as well as
a tsv file containg 0-based exclusive location of identical regions between 
the remaining genes of the pangenome. 
'''
# ---------------------------------------------------------------------------
from datetime import datetime
import os
import argparse
import subprocess
import pandas as pd
import logging
logging.basicConfig(level=logging.NOTSET)
import networkx as nx
import sys
sys.path.insert(1, '/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *

# ---------------------------------------------------------------------------
def runBlast(fastaPath):
	'''Run blast of the pangenome against itself with only exact matches (with minimum length of 100 bp)'''
	command = ['blastn', 
		'-query', fastaPath, 
		'-subject', fastaPath, 
		'-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen', 
		'-dust', 'no', 
		'-strand', 'plus', 
		'-word_size', '100', 
		'-penalty', '-10000', 
		'-ungapped', 
		'-out', 'PangenomeVSPangenome.Identical.blastn']
	subprocess.call(command)

def timenow():
	return(str(datetime.now()).split(' ')[1].split('.')[0])

# ---------------------------------------------------------------------------

# Initiate the parser
parser = argparse.ArgumentParser(description = '''
This script uses remove some genes of the pangenome that show high level of 
exact redundancy with other genes of the pangenome. 
Exact redundancy is find using blastn. Gene are removed in an iterative 
manner. The script output processed nucleotide and protein fastas, as well as
a tsv file containg 0-based exclusive location of identical regions between 
the remaining genes of the pangenome. 
'''
)

# Pangenome fasta (nucleotide)
parser.add_argument("-pn", "--PangenomeNucleotide", help="Fasta of the pangenome (nucl)", required=True)
# Pangenome fasta (protein)
parser.add_argument("-pp", "--PangenomeProtein", help="Fasta of the pangenome (prot)", required=True)

# Read arguments from the command line
args = parser.parse_args()
pangenomeNucleotidePath = os.path.abspath(args.PangenomeNucleotide)
outNucFasta = ".".join(pangenomeNucleotidePath.split(".")[0:-2])
pangenomeProteinPath = os.path.abspath(args.PangenomeProtein)
outProtFasta = ".".join(pangenomeProteinPath.split(".")[0:-2])

# Read the fasta file of the nucleotide pangenome
fastaPath = pangenomeNucleotidePath
fasta = Fasta(fastaPath)

### Find identical matches across the pangenome
logging.info(f'{timenow()}\tRunning blast of pangenome against itself')
runBlast(fastaPath)

### Process raw data
logging.info(f'{timenow()}\tProcess blast results')
# Read blast results
data = pd.read_csv("PangenomeVSPangenome.Identical.blastn", sep = "\t", names = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".split())
# Remove blast against itself
data = data.loc[data['qseqid'] != data['sseqid'],:]
# Add BED object for alignments
data.loc[:,'qBED'] = [BED(BEDcoordinates(data.loc[i,'qseqid'], data.loc[i,['qstart', 'qend']].min()-1, data.loc[i,['qstart', 'qend']].max())) for i in data.index]
# Add BED object for query gene
data.loc[:,'BED'] = [BED(BEDcoordinates(data.loc[i,'qseqid'], 0, data.loc[i,'qlen'])) for i in data.index]
rawData = data.copy()

# Start loop to remove genes with high level of redundancy iteratively
logging.info(f'{timenow()}\tStarting iteration to remove redundant genes')
nbGenesToRemove = 1
genesToRemove = []
while nbGenesToRemove > 0:

	### Remove genes to remove
	data = rawData.copy()
	data = data.loc[(~ data['qseqid'].isin(genesToRemove)) & (~ data['sseqid'].isin(genesToRemove)),:]

	### Cluster genes that have a shared sequence together
	G = nx.Graph()
	G.add_edges_from([(data.loc[index,'qseqid'], data.loc[index,'sseqid']) for index in data.index])
	# Get connected components
	clusters = list(nx.connected_components(G))

	### Compute coverage on each gene
	dataPerGene = data.groupby('qseqid').agg({'qBED': 'sum', 'qlen': 'first', 'BED': 'first'}).reset_index()
	dataPerGene.loc[:,'Cov'] = [dataPerGene.loc[i,'BED'].overlapLen(dataPerGene.loc[i,'qBED'], percent = False) for i in dataPerGene.index]
	dataPerGene.loc[:,'Specific'] = dataPerGene['qlen'] - dataPerGene['Cov'] # Length of the region that is specific to the gene, that is not shared with any other gene
	
	### Select genes to remove (having less than 100 specific bp, that is not shared with any other gene)
	allGenesToRemove = dataPerGene.loc[dataPerGene['Specific'] < 100,'qseqid'].tolist()
	nbGenesToRemove = len(allGenesToRemove)

	### Remove a single gene by connected component of the graph (because removing 1 gene could change the graph)
	genesToRemoveOnThisRun = []
	for c in clusters:
		genes = [x for x in allGenesToRemove if x in c]
		# If a gene has to be removed in the connected component, remove the last one in alphabetical order, in order to remove accessory genes first
		if len(genes) > 0:
			genesToRemoveOnThisRun += [sorted(genes)[-1]]
	genesToRemove += genesToRemoveOnThisRun
	logging.info(f'{timenow()}\t{nbGenesToRemove} gene(s) to remove, processing {len(genesToRemoveOnThisRun)} on this step')


logging.info(f'{timenow()}\tProcessing final fasta files. ')

# Get Nucleotide fasta file
fastaNuc = Fasta(pangenomeNucleotidePath)
fastaNuc2 = Fasta([seq for seq in fastaNuc.sequences if not seq.id in genesToRemove])
fastaNuc2.toFile(f'{outNucFasta}.NoRedundancy.cds.fna')

# Get Pep fasta file
fastaPep = Fasta(pangenomeProteinPath)
fastaPep2 = Fasta([seq for seq in fastaPep.sequences if not seq.id in genesToRemove])
fastaPep2.toFile(f'{outProtFasta}.NoRedundancy.pep.faa')

logging.info(f'{timenow()}\tExtract coordinates of identical sequences across the pangenome. ')

# Get coordinates of identical sequences in the nucleotide pangenome
runBlast(f'{outNucFasta}.NoRedundancy.cds.fna')
data = pd.read_csv("PangenomeVSPangenome.Identical.blastn", sep = "\t", names = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen".split())
data = data.loc[data['qseqid'] != data['sseqid'],:] # Remove blast against itself
data['qstart'] = data['qstart'] - 1 # Change to 0-based exclusive coordinates
data['sstart'] = data['sstart'] - 1 # Change to 0-based exclusive coordinates
data = data.loc[:,['qseqid', 'qstart', 'qend', 'sseqid', 'sstart', 'send']]
data.columns = ['Gene1', 'Start1', 'End1', 'Gene2', 'Start2', 'End2']
data.to_csv(f'{outNucFasta}.NoRedundancy.cds.IdenticalRegions.tsv', sep = "\t", index = False)

# Clean temporary files
logging.info(f'{timenow()}\tCleaning temporary files. ')
os.remove("PangenomeVSPangenome.Identical.blastn")

