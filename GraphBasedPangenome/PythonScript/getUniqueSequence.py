#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/04/27
# version ='4.0'
# ---------------------------------------------------------------------------
'''
This scripts takes as input the segments of a graph (fasta format), the list 
of segments originating from the reference genome, the blast results of the
segments on themselves, and outputs a fasta file of non-redundant segments, 
annotated as Ref or NonRef. Criteria to consider 2 segments as redundant is 
50% coverage and 95% sequence identity (option must be given when running the
blastn command). 
'''
# ---------------------------------------------------------------------------
import os
import argparse
import pandas as pd
import re
import sys
sys.path.insert(1, '/ccc/work/cont007/fg0006/loeglerv/GenomeAssemblyTools')
from Tools import *
import networkx as nx
import logging
logging.basicConfig(level=logging.INFO)
# ---------------------------------------------------------------------------
# Definitions

# ---------------------------------------------------------------------------

# =============
# GET ARGUMENTS
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = """
This scripts takes as input the segments of a graph (fasta format), the list 
of segments originating from the reference genome, the blast results of the
segments on themselves, and outputs a fasta file of non-redundant segments, 
annotated as Ref or NonRef. Criteria to consider 2 segments as redundant is 
50% coverage and 95% sequence identity (option must be given when running the
blastn command). 
""")

# Segments
parser.add_argument("-f", "--segmentsFasta", help="Segments of the graph at the fasta format", required=True)
# Ref segments
parser.add_argument("-r", "--refSegments", help="ID of the segments originating from the reference genome (one ID per line)", required=True)
# Self blast
parser.add_argument("-b", "--selfBlast", help="blastn of the segments on themselves, using -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'", required=True)

# Read arguments from the command line
args = parser.parse_args()
segmentFastaPath = os.path.abspath(args.segmentsFasta)
refSegmentsPath = os.path.abspath(args.refSegments)
selfBlastPath = os.path.abspath(args.selfBlast)

#segmentFastaPath = 'SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.fasta'
#refSegmentsPath = 'RefSegments.txt'
#selfBlastPath = 'SacePangenomeGraph.500Haplotypes.gfa.min0.1kb.Self.blastn'

logging.info('Importing data')
# Import segments fasta
segments = Fasta(segmentFastaPath)
# Import ref segments
refSegments = pd.read_csv(refSegmentsPath, names = ['ID']).iloc[:,0].tolist()
# Import blast of segments on themselves
data = pd.read_csv(selfBlastPath, sep = '\t', names = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs'.split(' '))
data = data.loc[data['qseqid'] != data['sseqid'],:]
data = data.loc[data['qcovs'] >= 50,:]

# Create segment redundancy graph
logging.info('Create sequence similarity graph')
G = nx.Graph()
nodes = segments.getID()
G.add_nodes_from(nodes)
edges = data.loc[:,['qseqid', 'sseqid']].drop_duplicates()
G.add_edges_from(list(zip(edges['qseqid'], edges['sseqid'])))

# Get the segment with the higher connectivity per connected component, then ordered by length
logging.info('Select the longuest segment per connected component')
lenDF = pd.DataFrame(data = {'Segment': segments.getID(), 'Len': segments.getLengths(), 'Sequences': segments.sequences}).sort_values(by = ['Len'], ascending = False)
representativeSegments = []
annotatedSegments = []
for c in nx.connected_components(G):
	clusterLenDF = lenDF.loc[lenDF['Segment'].isin(c),:].copy()
	if any(clusterLenDF['Segment'].isin(refSegments)):
		clusterLenDF['ID'] = clusterLenDF.loc[:,'Segment'] + '_InRef'
		clusterLenDF['NewSequences'] = [Sequence(clusterLenDF.loc[i,'ID'], clusterLenDF.loc[i,'Sequences'].seq) for i in clusterLenDF.index]
	else :
		clusterLenDF['ID'] = clusterLenDF['Segment'] + '_OutRef'
		clusterLenDF['NewSequences'] = [Sequence(clusterLenDF.loc[i,'ID'], clusterLenDF.loc[i,'Sequences'].seq) for i in clusterLenDF.index]
	representativeSegments += [clusterLenDF.iloc[0,4]]
	annotatedSegments += clusterLenDF['NewSequences'].tolist()

representativeSegments = Fasta(representativeSegments)
annotatedSegments = Fasta(annotatedSegments)

# Write fasta to file
representativeSegments.toFile(f'{re.sub(".fasta", "", segmentFastaPath)}.NoRedundancy.fasta')
annotatedSegments.toFile(f'{re.sub(".fasta", "", segmentFastaPath)}.Annotated.fasta')
logging.info(f'Output fasta written to {re.sub(".fasta", "", segmentFastaPath)}.NoRedundancy.fasta and {re.sub(".fasta", "", segmentFastaPath)}.Annotated.fasta')

