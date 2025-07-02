#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/03/06
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script filters the results of the blastp of non reference genes of Sace
pangenome on as custom database composed by RefSeq and other databases (Shen 
2018 - 332 Fungi species, Yue 2017 - 5 Sapa strains). It infers phylogenetic
origin of each gene based on sequence similarity. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import os
import csv
import argparse
import sys
sys.path.insert(1, '/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts/GenomeAssemblyTools')
#sys.path.insert(1, '/ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools')
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

# Read arguments from the command line
args = parser.parse_args()

pangenomeFastaPath=os.path.abspath(args.pangenomeFasta)
blastResultsPath=os.path.abspath(args.blastResults)
output=os.path.abspath(args.output)

#pangenomeFastaPath="/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa"
#blastResultsPath="/shared/home/vloegler/Pangenome_Sace1000ONT/04-analysis/08_BlastOnRefSeq/PangenomeNonRef_on_RefSeq_Shen2018_Yue2017.blastp"

# Get sequences IDs for the non ref genes of the pangenome
#pangenomeFastaPath='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains/Pangenome.NoRedundancy.pep.faa'
nonRefGenes = [x for x in Fasta(pangenomeFastaPath).getID() if x.startswith('YX')]

# Get taxid files for each categories
taxidfilesTuple = [('sace.txids', 4932), # Sace
		   ('sapa.txids', 27291), # Sapa
		   ('saccharomyces.txids', 4930), # Saccharomyces
		   ('sami.txids', 114525), # Sensu stricto
		   ('saku.txids', 114524), # Sensu stricto
		   ('saar.txids', 1160507), # Sensu stricto
		   ('saeu.txids', 1080349), # Sensu stricto
		   ('sauv.txids', 230603), # Sensu stricto
		   ('fungi.txids', 4751), # Fungi
		   ('eukaryote.txids', 2759), # Eukaryote
		   ('bacteria.txids', 2), # Bacteria
		   ('archaea.txids', 2157), # Archaea
		   ('viruses.txids', 10239)] # Viruses
for t in taxidfilesTuple:
	if not os. path. exists(t[0]): 
		os.system(f'get_species_taxids.sh -t {t[1]} > {t[0]}')

if not os. path. exists('sensustricto.txids'): 
	os.system('cat sace.txids sapa.txids sami.txids saku.txids saar.txids saeu.txids sauv.txids > sensustricto.txids')

sace = pd.read_csv("sace.txids", header = None).iloc[:,0].tolist()
sapa = pd.read_csv("sapa.txids", header = None).iloc[:,0].tolist()
sensustricto = pd.read_csv("sensustricto.txids", header = None).iloc[:,0].tolist()
saccharomyces = pd.read_csv("saccharomyces.txids", header = None).iloc[:,0].tolist()
fungi = pd.read_csv("fungi.txids", header = None).iloc[:,0].tolist()
eukaryote = pd.read_csv("eukaryote.txids", header = None).iloc[:,0].tolist()
bacteria = pd.read_csv("bacteria.txids", header = None).iloc[:,0].tolist()
archaea = pd.read_csv("archaea.txids", header = None).iloc[:,0].tolist()
viruses = pd.read_csv("viruses.txids", header = None).iloc[:,0].tolist()

blastResults = pd.read_csv(blastResultsPath, sep ='\t', header = None, names='qseqid sseqid stitle pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs sscinames scomnames staxids'.split())

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

# Find origin of each gene
# ========================
# Filter out Stenotrophomonas maltophilia because of obsolete records
#blastResults = blastResults.loc[blastResults['sscinames'] != 'Stenotrophomonas maltophilia',:]

# Take best hit for each query sequence
blastResults = blastResults.groupby(['qseqid']).first().reset_index()

# Add origin
blastResults.loc[:,'staxids'] = blastResults.loc[:,'staxids'].apply(lambda x: [int(y) for y in x.split(';')])
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in viruses for y in x]))), 'origin'] = "Virus"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in bacteria for y in x]))), 'origin'] = "Bacteria"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in archaea for y in x]))), 'origin'] = "Archaea"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in eukaryote for y in x]))), 'origin'] = "Eukaryote"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in fungi for y in x]))), 'origin'] = "Fungi"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in saccharomyces for y in x]))), 'origin'] = "Saccharomyces"
blastResults.loc[(blastResults['pident'] >= 93) & (blastResults['staxids'].apply(lambda x: any([y in sapa for y in x]))), 'origin'] = "Saccharomyces paradoxus"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in sace for y in x]))), 'origin'] = "Saccharomyces cerevisiae"
blastResults.loc[~(blastResults['staxids'].apply(lambda x: any([y in eukaryote+bacteria+archaea+viruses for y in x]))), 'origin'] = "Error"

# Add mechanism
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in eukaryote+bacteria+archaea+viruses for y in x]))), 'mechanism'] = "Candidate HGT"
blastResults.loc[(blastResults['pident'] >= 80) & (blastResults['staxids'].apply(lambda x: any([y in eukaryote for y in x]))), 'mechanism'] = "HGT"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in saccharomyces for y in x]))), 'mechanism'] = "Candidate introgression"
blastResults.loc[(blastResults['pident'] >= 80) & (blastResults['staxids'].apply(lambda x: any([y in saccharomyces for y in x]))), 'mechanism'] = "Introgression"
blastResults.loc[(blastResults['staxids'].apply(lambda x: any([y in sace for y in x]))), 'mechanism'] = "Fast evolving gene"
blastResults.loc[~(blastResults['staxids'].apply(lambda x: any([y in eukaryote+bacteria+archaea+viruses for y in x]))), 'mechanism'] = "Error"

# Add genes with no sequence similarity in the Database
noHit = pd.DataFrame(columns = blastResults.columns).astype(blastResults.dtypes)
noHit.loc[:,'qseqid'] = [x for x in nonRefGenes if not x in blastResults['qseqid'].tolist()]
noHit.loc[:,'origin'] = 'NoHit'
noHit.loc[:,'mechanism'] = 'Unknown'

origins = pd.concat([blastResults, noHit])
origins.to_csv(output, sep = '\t', index = False, quoting = csv.QUOTE_NONNUMERIC)
