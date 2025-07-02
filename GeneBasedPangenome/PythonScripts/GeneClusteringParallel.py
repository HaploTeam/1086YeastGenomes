#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/04/27
# version ='4.0'
# ---------------------------------------------------------------------------
'''
This script uses a graph-based approach to create groups of orthologous genes
from denovo genome annotations. The clustering is based on All VS All blastn
of non reference genes (redundant) plus reference genes. 

Here is the blast command used for each strain:
	blastn  -query $WORKDIR/03_NonOrthoGenes/$S.NonOrthoGenes.cds.fa 
			-db AllNonOrthoGenes.fa 
			-num_threads 4 
			-perc_identity 95
			-strand plus 
			-dust no
			-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' 
			-out $S.vsAll.blastn

Once the genes clustering is done, a representative gene is chosen for each 
groups of orthologs and fasta files are written. 
'''
# ---------------------------------------------------------------------------
import re
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1' # Required to avoid "OpenBLAS blas_thread_init: pthread_create: Resource temporarily unavailable" error
import time
from datetime import datetime
import argparse
import pandas as pd
import sys
sys.path.insert(1, '/shared/home/vloegler/Pangenome_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import networkx.algorithms.community as nx_comm
import logging
logging.basicConfig(level=logging.INFO)
import itertools
from multiprocessing import Pool
import pickle
# ---------------------------------------------------------------------------
# Definitions
def timenow():
	return(str(datetime.now()).split(' ')[1].split('.')[0])

# ---------------------------------------------------------------------------

# =============
# GET ARGUMENTS
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = """
This script uses a graph-based approach to create groups of orthologous genes
from denovo genome annotations. The clustering is based on All VS All blastn
of non reference genes (redundant) plus reference genes. 

Here is the blast command used for each strain:
	blastn  -query $WORKDIR/03_NonOrthoGenes/$S.NonOrthoGenes.cds.fa 
			-db AllNonOrthoGenes.fa 
			-num_threads 4 
			-perc_identity 95 
			-strand plus 
			-dust no
			-outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen' 
			-out $S.vsAll.blastn

Once the genes clustering is done, a representative gene is chosen for each 
groups of orthologs and fasta files are written. 
""")

# All genes fasta
parser.add_argument("-a", "--allGenesFasta", help="Fasta of all genes with no orthology to the reference (nucl)", required=True)
# Reference genes
parser.add_argument("-rg", "--referenceGenesFasta", help="Fasta of all genes of the reference", required=True)
parser.add_argument("-rp", "--referenceProteinsFasta", help="Fasta of all proteins of the reference", required=True)
# Transfered annotations directory
parser.add_argument("-annot", "--annotationDir", help="Directory where are stored the annotation, fasta cds, trimmed_cds and pep for each strain", required=True)
# Directory containing ref genes present in each strain
parser.add_argument("-ortho", "--orthoGenesDir", help="Directory where are stored the fasta files of the genes with an orthology to the reference", required=True)
# Directory containing results of Blast All VS All non ref genes
parser.add_argument("-b", "--blastAllVSAllDir", help="Directory containing the blastn results of All vs All genes", required=True)
# Output directory
parser.add_argument("-o", "--outputDir", help="Output directory", required=True)
# Number of threads
parser.add_argument("-t", "--nbThreads", help="Number of threads for parallelization", default = 1, type = int)
# Output localization dataframe (dataframe containing localization information about non-reference genes) and Fasta sequences in addition to pangenome matrix
parser.add_argument('-n', '--noCompleteOutput', help="Do not output localization dataframe (dataframe containing localization information about non-reference genes) and Fasta sequences in addition to pangenome matrix", action='store_true')

# Read arguments from the command line
args = parser.parse_args()
nbThreads = args.nbThreads
noCompleteOutput = args.noCompleteOutput
allGenesFastaPath = os.path.abspath(args.allGenesFasta)
referenceGenesFastaPath = os.path.abspath(args.referenceGenesFasta)
referenceProteinsFastaPath = os.path.abspath(args.referenceProteinsFasta)
annotationDir = os.path.abspath(args.annotationDir)
if annotationDir.endswith('/'):
	annotationDir = annotationDir[:-1]
orthoGenesDir = os.path.abspath(args.orthoGenesDir)
if orthoGenesDir.endswith('/'):
	orthoGenesDir = orthoGenesDir[:-1]
blastAllVSAllDir = os.path.abspath(args.blastAllVSAllDir)
if blastAllVSAllDir.endswith('/'):
	blastAllVSAllDir = blastAllVSAllDir[:-1]
outputDir = os.path.abspath(args.outputDir)
if outputDir.endswith('/'):
	outputDir = outputDir[:-1]

# Create output directory if it does not exists, and GraphPlots subdirectory
if not os.path.exists(outputDir):
	logging.info(f'{timenow()}\tCreate output directory')
	os.makedirs(outputDir)
if not os.path.exists(f'{outputDir}/GraphPlots'):
	os.makedirs(f'{outputDir}/GraphPlots')

#allGenesFastaPath='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/05_BlastNonOrthoAllVSAll_HQStrains/AllGenes_Redundant.min0.1kb.fasta'
#referenceGenesFastaPath='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/05_BlastNonOrthoAllVSAll_HQStrains/orf_genomic_all_R64-3-1_20210421.min0.1kb.fasta'
#referenceProteinsFastaPath='/home/vloegler/Pangenome_Sace1000ONT/03-data/orf_trans_all_R64-3-1_20230313.fasta'
#annotationDir='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/01_TransferAnnotation'
#orthoGenesDir='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/02_RefGenes'
#blastAllVSAllDir='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/05_BlastNonOrthoAllVSAll_HQStrains'
#outputDir='/home/vloegler/Pangenome_Sace1000ONT/04-analysis/07_GeneClustering_HQStrains'
#nbThreads=4

# Get ID of all sequences
logging.info(f'{timenow()}\tRead gene fasta files')
completeFasta = Fasta(allGenesFastaPath)
referenceFasta = Fasta(referenceGenesFastaPath)
IDs = pd.Series(completeFasta.getID())
nonRefIDs = IDs.loc[~IDs.isin(referenceFasta.getID())]

# Get all haplotypes from the list of blastn results files
haplotypes = os.listdir(blastAllVSAllDir)
haplotypes = sorted([re.sub('.vsAll.blastn', '', x) for x in haplotypes if x.endswith('.vsAll.blastn') and x != 'Reference.vsAll.blastn'])

# ================
# NETWORK BUILDING
# ================

# Create a Networkx graph
# Each node is a protein
# There is an edge between 2 nodes if a blast hit match pass the thresholds
G = nx.Graph()

def checkOrthology(blastResultsPath):
	'''This function outputs pairs of genes that are orthologous by
	analysing blast results'''
	# Inspect blast results
	blastOut = pd.read_csv(blastResultsPath, sep = '\t', header = None, names = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen'.split())
	# Remove lines where subject and query are the same gene
	blastOut = blastOut.loc[blastOut['qseqid'] != blastOut['sseqid'],:]
	# Convert Alignments and genes into BED objects
	blastOut.loc[:,'sBED'] = blastOut['slen'].apply(lambda x: BED(BEDcoordinates('Gene', 0, x)))
	blastOut.loc[:,'qBED'] = blastOut['qlen'].apply(lambda x: BED(BEDcoordinates('Gene', 0, x)))
	blastOut.loc[:,'sA_BED'] = blastOut.loc[:,['sstart', 'send']].apply(lambda x: BED(BEDcoordinates('Gene', x.iloc[0]-1, x.iloc[1])), axis = 1)
	blastOut.loc[:,'qA_BED'] = blastOut.loc[:,['qstart', 'qend']].apply(lambda x: BED(BEDcoordinates('Gene', x.iloc[0]-1, x.iloc[1])), axis = 1)
	# Group by QuerySubject pairs and sum alignment BEDs
	blastOutGrouped = blastOut.groupby(['qseqid', 'sseqid']).agg({'sBED': 'first', 'qBED': 'first', 'sA_BED':'sum','qA_BED':'sum'}).reset_index()
	# Compute query and subject coverage
	blastOutGrouped.loc[:,'scov'] = blastOutGrouped.apply(lambda x: x.sBED.overlapLen(x.sA_BED, percent = True), axis = 1)
	blastOutGrouped.loc[:,'qcov'] = blastOutGrouped.apply(lambda x: x.qBED.overlapLen(x.qA_BED, percent = True), axis = 1)
	# Keep match with query and subject coverage superior to 50, or one coverage superior to 90
	blastOutGrouped = blastOutGrouped.loc[
		((blastOutGrouped['qcov'] >= 50) & (blastOutGrouped['scov'] >= 50)) | 
		(blastOutGrouped['scov'] >= 90) |
		(blastOutGrouped['qcov'] >= 90),
	:]
	# Return list of tuples corresponding to orthologous genes
	return([(blastOutGrouped.loc[index,'qseqid'], blastOutGrouped.loc[index,'sseqid']) for index in blastOutGrouped.index])

# Check orthology over each blastn results files
# ==============================================
# List blastn results files
blastResultsFiles = [f'{blastAllVSAllDir}/{x}' for x in os.listdir(blastAllVSAllDir) if x.endswith('.vsAll.blastn')]

# Run orthology check in parallel
# Create the pool
process_pool = Pool(processes=nbThreads)
# Start processes in the pool
logging.info(f'{timenow()}\tCheck for orthology relationships in {len(blastResultsFiles)} blast result files')
start_time = time.time()
orthologyList = process_pool.map(checkOrthology, blastResultsFiles)
process_pool.close()
logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')
# Merge list of lists into simple list
orthologyList = list(itertools.chain.from_iterable(orthologyList))

# Add edges to the network
# ========================
G.add_edges_from(orthologyList)

# Write the network to file
# =========================
nx.write_multiline_adjlist(G,f"{outputDir}/GeneClusteringGraph.multiline-adjlist")

# To skip first par of the script, just read the network file
#G = nx.read_multiline_adjlist(f"{outputDir}/GeneClusteringGraph.multiline-adjlist")


# ===============
# GENE CLUSTERING
# ===============

# 1. Defining clusters ========================================================

# Each connected component of the graph is considered as a gene cluster. 
# If the component density is lower than 0.4, Louvain's communities are 
# detected and each community is considered as a cluster. When density
# is lower than 0.4, the component is plotted with different colors for 
# the different communities

# 2. Finding the representative gene of each cluster ==========================

# The best representative gene is the node with the higher degree (number of
# edges). If mutiple nodes have the same degree, the longuest gene is taken. If
# a reference gene is present in the orthology group, it is taken as
# representative. 

def defineOrthologyGroup(cluster):
	'''
	This function takes a list of gene that correponds to a connected component 
	of the graph (cluster) and defines the orthology groups within this cluster.
	- If the density of the cluster in the graph is higher 
	than 0.4, the cluster corresponds to a single orthology group and the 
	representative gene for the group is defined as the most connected and 
	longer reference gene, or non reference gene if no reference gene is 
	present. 
	- When the density is lower than 0.4, the cluster is further split into 
	Louvain communities, which all corresponds to a group of orthology. The 
	representative gene for each group is chosen as described above. 

	The function output a list of tuples: 
	[(RepresentativeGene, [list of genes in the orthology group]), 
	(RepresentativeGene, [list of genes in the orthology group]), 
	...]
	'''
	# Get a subraph with only genes present in the cluster
	subG = G.subgraph(cluster)
	# Get the cluster density
	d = nx.density(subG)
	# Define groups of orthology in the cluster
	# =========================================
	if d >= 0.4: # If density is greater or equal to 0.4, the cluster corresponds to a single group or orthology
		# Cluster is considered as a single community
		communities = [cluster]
	else: # If density is lower than 0.4
		# Find louvain communities
		communities = nx_comm.louvain_communities(subG, seed=123)
		# Plot the graph
		plt.figure(figsize=(8,8))
		plt.axis("off")
		plt.title(f'Connected component with {subG.number_of_nodes()} nodes, {subG.number_of_edges()} edges and {len(communities)} communities')
		options = {"edgecolors": "tab:gray", "node_size": 200, "alpha": 0.9}
		pos = nx.spring_layout(subG, seed=3113794652)  # positions for all nodes
		nx.draw_networkx_edges(subG, pos, width=1.0, alpha=0.5) # Draw edges
		for comm in communities: # Draw nodes of each Louvain communities with a random color
			rgb = np.random.rand(3,)
			nx.draw_networkx_nodes(subG, pos, nodelist=comm, node_color=[rgb], **options)
		plt.savefig(f'{outputDir}/GraphPlots/Cluster_{sorted(list(subG.nodes()))[0]}.png') # Save plot to png
		plt.close()
	# Find the representative gene of each community
	outputList = []
	for comm in communities:
		# Get subgraph for community
		subSubG = subG.subgraph(comm)
		# Check if a reference gene is present
		refGenesPresent = [x for x in referenceFasta.getID() if x in comm]
		if len(refGenesPresent) > 0:
			# Find the reference gene with the higher degree
			maxDegree = max([subSubG.degree(x) for x in refGenesPresent])
			nodesMaxDeg = [x for x in refGenesPresent if subSubG.degree(x) == maxDegree]
			if len(nodesMaxDeg) == 1: # If a single node with max degree
				representativeGene = nodesMaxDeg[0]
			else: # If multiple nodes with max degree, take the longuest gene
				# Retrieve gene length
				lengths = [len(referenceFasta.getSeqFromID(x)) for x in nodesMaxDeg]
				representativeGene = [nodesMaxDeg[j] for j in range(len(nodesMaxDeg)) if lengths[j] == max(lengths)][0]
		else: # If no reference gene is present
			# Find nodes with maximum degree
			maxDegree = max([subSubG.degree(x) for x in comm])
			nodesMaxDeg = [x for x in comm if subSubG.degree(x) == maxDegree]
			if len(nodesMaxDeg) == 1: # If a single node with max degree
				representativeGene = nodesMaxDeg[0]
			else: # If multiple nodes with max degree, take the longuest gene
				# Retrieve gene length
				lengths = [len(completeFasta.getSeqFromID(x)) for x in nodesMaxDeg]
				representativeGene = [nodesMaxDeg[j] for j in range(len(nodesMaxDeg)) if lengths[j] == max(lengths)][0]
		# Increment output list with a tuple containing the representative gene and the list of genes contained in the group of orthology
		outputList += [(representativeGene, list(subSubG.nodes()))]
	return(outputList)


# Define orthology groups from each connected component of the graph
# ==================================================================

# Separate each connected component of the graph
clusters = list(nx.connected_components(G))

# Run orthology group definition in parallel
# Create the pool
process_pool = Pool(processes=nbThreads)
# Start processes in the pool
logging.info(f'{timenow()}\tDefine orthology groups for {len(clusters)} connected components in the graph')
start_time = time.time()
orthologyGroups = process_pool.map(defineOrthologyGroup, clusters)
process_pool.close()
logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')
# Merge list of lists into simple list
orthologyGroups = list(itertools.chain.from_iterable(orthologyGroups))

# Create dictionary from list of tuples
orthologyGroups = dict(orthologyGroups) # Dictionary that will contain key(Representative gene) = [list of genes in cluster]

# Retrieve singleton genes (genes that didn't cluster with anything)
logging.info(f'{timenow()}\tAdd singleton genes to the graph')
start_time = time.time()
singleGenes = [x for x in IDs if not x in list(G.nodes())]
# Adding singleton genes (genes present in a single haplotype), and ref genes that didn't cluster with anything
for gene in singleGenes:
	orthologyGroups[gene] = [gene]
logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')

# ====================================
# ADD REFERENCE GENES ORTHOLOGY GROUPS
# ====================================

# Since only genes with no orthology to the reference genes were used for the 
# gene clustering, we need to build orthology groups for reference genes. We
# will use information from the reference gene annotation transfer. 

# Add information from reference gene annotation transfer
# =======================================================

logging.info(f'{timenow()}\tAdd orthology information from reference annotation transfer for {len(haplotypes)} haplotypes')
start_time = time.time()

# Import denovo genes with reference annotation transfer
for h in haplotypes:
	fasta = Fasta(f'{orthoGenesDir}/{h}.OrthoGenes.cds.fa')
	# Genes are stored in a list of tuples (ReferenceGeneName, HaplotypeGeneID)
	genes = [(x.split('|')[2], x.split('|')[0]) for x in fasta.getID()]
	# Remove truncated information
	genes = [(re.sub('_truncated', '', x[0]), x[1]) for x in genes]
	# Add gene presence to the pangenome dictionary
	for gene in genes:
		for g in gene[0].split('/'): # Split with '/' for cases where a gene have multiple orthologs in the reference
			if g in orthologyGroups.keys() and not gene[1] in orthologyGroups[g]:
				orthologyGroups[g] = list(orthologyGroups[g]) + [gene[1]]

logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')

# Write dictionary to file
with open(f'{outputDir}/OrthologyDictionary.pkl', 'wb') as fp:
	pickle.dump(orthologyGroups, fp)
	logging.info(f'{timenow()}\tOrthology dictionary saved to file')

# Load orhtology dictionary
#with open(f'{outputDir}/OrthologyDictionary.pkl', 'rb') as fp:
#	orthologyGroups = pickle.load(fp)


# ==========================
# CREATE PANGENOME DATAFRAME
# ==========================

# Get genes present in each group of orthology for each haplotype
# ===============================================================

### Convert dictionnary to dataframe with all pairs of RepresentativeGene, Gene
# Create list of tuples with (key, gene)
orthologyGroupsDF = [[(key, gene) for gene in orthologyGroups[key]] for key in orthologyGroups.keys()]
# Merge list of lists into simple list
orthologyGroupsDF = list(itertools.chain.from_iterable(orthologyGroupsDF))
# Convert to dataframe
orthologyGroupsDF = pd.DataFrame(orthologyGroupsDF, columns =['Key', 'Gene'])

def getGenePresence(haplotype):
	'''
	This function outputs a pandas DataFrame object that contains the gene 
	presence information for each orthology group in the haplotype. If one or 
	many genes are present in the haplotype for a given orthology group, value 
	will be 'Gene1,Gene2,Gene3'. If no gene is present in the haplotype for a 
	given orthology group, value will be '*'. 
	'''
	# Get the total list of genes present in the haplotype from the CDS fasta file of the haplotype
	genes = Fasta(f'{annotationDir}/{haplotype}.nuclear_genome.transAnnot.cds.fa').getID()
	genes = [x.split('|')[0] for x in genes] # Take only haplotype specific gene ID
	# Filter orthology gourps dataframe for the haplotype
	haploDF = orthologyGroupsDF.rename(columns = {'Key': 'Gene', 'Gene': haplotype})
	haploDF = haploDF.loc[haploDF[haplotype].isin(genes),:]
	# # Group genes belonging to the same orthology group
	haploDF = haploDF.groupby('Gene').agg(','.join)
	# Merge haplotype specific genes with all groups dataframe
	completeDF = pd.DataFrame(data = {'Gene': list(orthologyGroups.keys())})
	completeDF = completeDF.merge(haploDF, on = 'Gene', how = 'left')
	# Set gene to '*' when no gene is present
	completeDF.loc[completeDF[haplotype].isnull(), haplotype] = '*'
	# Set index and drop Gene column
	completeDF.index = completeDF['Gene']
	completeDF = completeDF.drop(columns = ['Gene'])
	return(completeDF)

# Get gene presence information for each haplotype in parallel
# Create the pool
process_pool = Pool(processes=nbThreads)
# Start processes in the pool
logging.info(f'{timenow()}\tBuild pangenome dataframe for {len(haplotypes)} haplotypes')
start_time = time.time()
pangenomePhased = process_pool.map(getGenePresence, haplotypes)
process_pool.close()
logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')
# Concatenate pangenomePhased in a single dataframe
pangenomePhased = pd.concat(pangenomePhased, axis = 1)

# Merge haplotypes together
# =========================

# Get strains from haplotypes
strains = sorted(set([re.sub('.HP1', '', x) for x in [re.sub('.HP2', '', y) for y in haplotypes]]))

def mergeHaplotypes(strain):
	'''
	This function takes as input a strain, merge haplotypes from the 
	pangenomePhased columns if the strain is phased and outputs a dataframe 
	with a single column corresponding to all the haplotypes of the strain. 
	'''
	# Get all different haplotypes
	hps = [x for x in haplotypes if x == f'{strain}.HP1' or x == f'{strain}.HP2' or x == strain]
	# Get pangenomePhased columns of these haplotypes
	df = pangenomePhased.loc[:,hps]
	# If multiple haplotypes present, merge them
	if len(hps) > 1:
		# Merge columns
		df[strain] = df[hps].agg(','.join, axis=1)
		df[strain] = df[strain].apply(lambda x: re.sub('\*,\*', '*', x))
		df[strain] = df[strain].apply(lambda x: re.sub(',\*', '', x))
		df[strain] = df[strain].apply(lambda x: re.sub('\*,', '', x))
		# Drop haplotype columns
		df = df.drop(columns = hps)
	return(df)

# Merge haplotypes for each strain in parallel
# Create the pool
process_pool = Pool(processes=nbThreads)
# Start processes in the pool
logging.info(f'{timenow()}\tMerge pangenome of {len(haplotypes)} haplotypes in {len(strains)} strains')
start_time = time.time()
pangenome = process_pool.map(mergeHaplotypes, strains)
process_pool.close()
logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')
# Concatenate pangenome in a single dataframe
pangenome = pd.concat(pangenome, axis = 1).copy()

# Add reference strain to the dataframe
# =====================================

logging.info(f'{timenow()}\tAdd the reference strain S288c to the pangenome')
start_time = time.time()

# Get the total list of genes present in the reference strain from the CDS fasta file of the strain
genes = referenceFasta.getID()
# Filter dataframe containing orthology relationship to keep only refrence genes
haploDF = orthologyGroupsDF.rename(columns = {'Key': 'Gene', 'Gene': 'S288c'})
haploDF = haploDF.loc[haploDF['S288c'].isin(genes),:]
# Group genes belonging to the same orthology group
haploDF = haploDF.groupby('Gene').agg(','.join)
# Merge reference genes with all groups dataframe
refPangenome = pd.DataFrame(data = {'Gene': list(orthologyGroups.keys())})
refPangenome = refPangenome.merge(haploDF, on = 'Gene', how = 'left')
# Set gene to '*' when no gene is present
refPangenome.loc[refPangenome['S288c'].isnull(), 'S288c'] = '*'
# Set index and drop Gene column
refPangenome.index = refPangenome['Gene']
refPangenome = refPangenome.drop(columns = ['Gene'])

# Merge to pangenome
pangenome = pd.concat([pangenome, refPangenome], axis = 1)
logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')

# Output pangenome to file
# ========================
pangenome.to_csv(f'{outputDir}/Pangenome.tsv', sep = '\t')
pangenomePhased.to_csv(f'{outputDir}/PangenomePhased.tsv', sep = '\t')
logging.info(f'{timenow()}\tPangenome written to {outputDir}/Pangenome.tsv')

# Write tsv file with reference gene families
refGenes = pangenome.copy()
refGenes['RepresentativeGene'] = refGenes.index
refGenes = refGenes.loc[:,['RepresentativeGene', 'S288c']]
refGenes = refGenes.loc[refGenes['S288c'] != '*',:].rename(columns = {'S288c': 'ReferenceGenes'})
refGenes.to_csv(f'{outputDir}/ReferenceGeneFamilies.tsv', sep = '\t', index = False)

# ========================
# LOCALIZATION INFORMATION
# ========================
# Find localization informations for non reference genes

if not noCompleteOutput:

	def getLocalizationInfo(ID):
		'''
		This function finds the coordinates of a given gene, as well as the 
		upstream and downstream genes.
		'''
		# Find pangenomePhased column containing ID
		subPangenome = pangenomePhased.loc[:, pangenomePhased.loc[ID,:].str.contains(ID)]
		# Deduct the strain containing the gene
		s = subPangenome.columns.tolist()[0]
		# Read gff3 file
		annot = pd.read_csv(f'{annotationDir}/{s}.nuclear_genome.transAnnot.gff3', sep = "\t", comment="#", header = None, names = "seqid source type start end score strand phase attributes".split())
		annot = annot.loc[annot['type'] == "gene",:]
		# Extract gene ID and name
		annot['ID'] = [x.split(';')[0].split('=')[1] for x in annot['attributes']]
		annot['name'] = [x.split('=')[2] for x in annot['attributes']]
		# Keep only genes that were not considered for the pangenome construction (ex: genes larger than 100 bp) or genes with orthology to the reference (meaning that ID != Name)
		annot = annot.loc[(annot['ID'] != annot['name']) | (annot['ID'].isin(completeFasta.getID())),:]
		# Get chromosome, strand, start and end information
		chromosome = annot.loc[annot['attributes'] == f'ID={ID};Name={ID}','seqid'].iloc[0]
		strand = annot.loc[annot['attributes'] == f'ID={ID};Name={ID}','strand'].iloc[0]
		start = annot.loc[annot['attributes'] == f'ID={ID};Name={ID}','start'].iloc[0]
		end = annot.loc[annot['attributes'] == f'ID={ID};Name={ID}','end'].iloc[0]
		# Get next and previous genes on the chromosome, None if the gene is first or last
		chrAnnot = annot.loc[annot['seqid'] == chromosome,:].reset_index()
		index = chrAnnot.loc[chrAnnot['attributes'] == f'ID={ID};Name={ID}',:].index
		if index == 0:
			previousGene = None
		else:
			previousGene = chrAnnot.iloc[index - 1,9].iloc[0].split("=")[2]
		if index == max(chrAnnot.index):
			nextGene = None
		else:
			nextGene = chrAnnot.iloc[index + 1,9].iloc[0].split("=")[2]
		return(pd.DataFrame({'RepresentantGeneID': [ID], 'OriginStrain': [s], 'Chromosome': [chromosome], 'Strand': [strand], 'Start': [start], 'End': [end], 'PreviousGene': [previousGene], 'NextGene': [nextGene]}))

	# Get localization information for each gene in parallel
	nonRefGenes = [x for x in pangenome.index if not x in referenceFasta.getID()]

	# Create the pool
	process_pool = Pool(processes=nbThreads)
	# Start processes in the pool
	logging.info(f'{timenow()}\tGet localization information for {len(nonRefGenes)} non reference genes')
	start_time = time.time()
	locDF = process_pool.map(getLocalizationInfo, nonRefGenes)
	process_pool.close()
	logging.info(f'{timenow()}\tFinished in {round(time.time() - start_time, 1)} seconds')
	# Concatenate pangenome in a single dataframe
	locDF = pd.concat(locDF, axis = 0).reset_index(drop = True)

	# Order by chromosome, then by start position
	def getChrNumber(chrID):
		chr = [x for x in chrID.split("_") if 'chromosome' in x or 'Unplaced' in x][0]
		if 'chromosome' in chr:
			return int(re.sub('chromosome', '', chr.split('.')[0]))
		else:
			return 17 # return a number higher than the last chromosome

	locDF.loc[:,'CHR'] = locDF['Chromosome'].apply(getChrNumber)
	locDF = locDF.sort_values(by = ['CHR', 'Start']).reset_index(drop = True)

	# Add standardized geneID
	newNames = dict([(locDF.iloc[i,locDF.columns.get_loc('RepresentantGeneID')], f'YX{(1 + i):05d}') for i in locDF.index])
	locDF.loc[:,'GeneID'] = locDF['RepresentantGeneID'].apply(lambda x: newNames[x])
	locDF = locDF.loc[:,['GeneID', 'RepresentantGeneID', 'OriginStrain', 'Chromosome', 'Strand', 'Start', 'End', 'PreviousGene', 'NextGene']]

	# Change ID in pangenome dataframe
	oldIndex = pangenome.index
	newIndex = [newNames[x] if x in newNames.keys() else x for x in pangenome.index]
	pangenome.index = newIndex
	# Order by lexicographical order
	pangenome = pangenome.loc[sorted(pangenome.index),:]


	# Get a dataframe containing all genes for each orthology group in a single column
	pangenomeSingleCol = pangenome.agg(','.join, axis = 1)

	# Change previous and next gene to standardized ID
	def changeID(x):
		# if gene is not a reference gene
		if x != None and not x in referenceFasta.getID() and x in completeFasta.getID():
			# Find pangenome row containing ID
			x = pangenomeSingleCol[pangenomeSingleCol.str.contains(x)].index.values[0]
		return(x)

	locDF.loc[:,'PreviousGene'] = locDF['PreviousGene'].apply(changeID)
	locDF.loc[:,'NextGene'] = locDF['NextGene'].apply(changeID)
	# Write to file ====================================================================
	pangenome.to_csv(f'{outputDir}/Pangenome.tsv', sep = '\t')
	locDF.to_csv(f'{outputDir}/Pangenome.NonRef.Localization.tsv', sep = "\t", index = False)
	logging.info(f'{timenow()}\tLocalization information written to {outputDir}/Pangenome.NonRef.Localization.tsv')

	# ========================
	# RETRIEVE FASTA SEQUENCES
	# ========================

	# Create pangenome nucleotide fasta file
	# ======================================
	## Convert Fasta object to dataframe to filter quickly sequences present in pangenome.index
	dfCompleteFasta = pd.DataFrame({'SeqObjects': completeFasta.sequences})
	dfCompleteFasta.index = dfCompleteFasta['SeqObjects'].apply(lambda x: x.id)
	dfPangenomeFasta = dfCompleteFasta.loc[oldIndex,:]
	# Change ID of the non ref genes with the new name
	def changeID(seqObject, newID):
		seqObject.id = newID
		seqObject.description = f'>{newID}\n'
		return(seqObject)
	dfPangenomeFasta.loc[:,'SeqObjects'] = [changeID(seq, newNames[seq.id]) if seq.id in newNames.keys() else seq for seq in dfPangenomeFasta.loc[:,'SeqObjects']]
	# Order by lexicographical order
	dfPangenomeFasta.index = dfPangenomeFasta['SeqObjects'].apply(lambda x: x.id)
	dfPangenomeFasta = dfPangenomeFasta.loc[sorted(dfPangenomeFasta.index),:]
	# Write to file
	pangenomeFasta = Fasta(dfPangenomeFasta['SeqObjects'].tolist())
	pangenomeFasta.toFile(f'{outputDir}/Pangenome.cds.fna')
	logging.info(f'{timenow()}\tPangenome cds written to {outputDir}/Pangenome.cds.fna')

	# Create pep fasta file
	# =====================
	## Get sequences of reference strain, then sequences for each strains
	## Then convert to fasta

	### Add reference sequences
	# Read reference pep fasta
	referencePepFasta = Fasta(referenceProteinsFastaPath)
	# Retrieve sequences that are representative of an orthology group
	pepSequences = [seq for seq in referencePepFasta.sequences if seq.id in oldIndex]

	### Add sequences for each strains
	def getRepresentativeGenePepSeq(haplotype):
		'''
		This function retrieve all the peptide sequences corresponding to a 
		representative gene of the pangenome for a given strain
		'''
		# Read peptide fasta for the strain
		pepFasta = Fasta(f'{annotationDir}/{haplotype}.nuclear_genome.transAnnot.pep.fa')
		# Change IDs and description of sequences
		for seq in pepFasta.sequences:
			seq.id = seq.id.split('|')[0]
			seq.description = f'>{seq.id}\n'
		# Get genes that are representative of an orthology group
		genes = [x for x in pepFasta.getID() if x in oldIndex]
		# If some genes are found, output the corresponding sequence objects
		if len(genes) > 0:
			# Retrieve seq objects
			seq = [s for s in pepFasta.sequences if s.id in genes]
			return(seq)
		else:
			return([])

	for haplotype in haplotypes:
		pepSequences += getRepresentativeGenePepSeq(haplotype)

	# Change ID of sequence object that are not reference proteins
	pepSequences = [changeID(seq, newNames[seq.id]) if seq.id in newNames.keys() else seq for seq in pepSequences]

	# Order by lexicographical order
	dfPepFasta = pd.DataFrame({'SeqObjects': pepSequences})
	dfPepFasta.index = dfPepFasta['SeqObjects'].apply(lambda x: x.id)
	dfPepFasta = dfPepFasta.loc[sorted(dfPepFasta.index),:]
	# Write to file
	pangenomePepFasta = Fasta(dfPepFasta['SeqObjects'].tolist())
	pangenomePepFasta.toFile(f'{outputDir}/Pangenome.pep.faa')
	logging.info(f'{timenow()}\tPangenome pep sequences written to {outputDir}/Pangenome.cds.fna')

