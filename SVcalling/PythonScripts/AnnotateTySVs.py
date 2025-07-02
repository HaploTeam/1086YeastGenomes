#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/10/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script takes a fasta of SV sequences and detects Ty-related SVs.
'''
# ---------------------------------------------------------------------------
import subprocess
import pandas as pd
import argparse
import os
import logging
import subprocess
import sys
sys.path.insert(1, '/home/vloegler/SVCalling_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
logging.basicConfig(level=logging.NOTSET)
# ---------------------------------------------------------------------------
# Definitions
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--fasta", help="Path of the SV sequences (fasta)", required=True)
parser.add_argument("-ty", "--ty", help="Path of the Ty sequences (fasta)", required=True)
parser.add_argument("-o", "--output", help="Output file (VCF format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

fastaPath=os.path.abspath(args.fasta)
tyPath=os.path.abspath(args.ty)
outputPath=os.path.abspath(args.output)

#fastaPath='SV.1087Samples.DeletionSequences.fasta'
#tyPath='TyDB.fasta'

# Blast Ty on SVs ============================================================

# Make blast db
cmd = ['makeblastdb', 
       '-dbtype', 'nucl', 
       '-title', 'TyDB', 
       '-out', 'TyDB', 
       '-in', 'TyDB.fasta']
subprocess.call(cmd)

# Run blast against Ty DB
cmd = ['blastn', 
       '-query', fastaPath, 
       '-db', 'TyDB', 
       '-dust', 'no', 
       '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen', 
       '-out', 'TySearch.blastn']
subprocess.call(cmd)

# Process blast results =====================================================
results = pd.read_csv('TySearch.blastn', sep = '\t', names = 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen'.split(' '))

# Filter hsp with less than 95% identity
results = results.loc[results['pident'] >= 95,:]

# Add BED object of the hsp
results['alnBED'] = [BED(BEDcoordinates(results.loc[i,'qseqid'], results.loc[i,['qstart', 'qend']].min() - 1, results.loc[i,['qstart', 'qend']].max())) for i in results.index]

# Group by SV and add BED objects
resultsGroupped = results.groupby(['qseqid', 'qlen'])['alnBED'].sum().reset_index()

# Add BED object of the SV
resultsGroupped['svBED'] = [BED(BEDcoordinates(resultsGroupped.loc[i,'qseqid'], 0, resultsGroupped.loc[i,'qlen'])) for i in resultsGroupped.index]

# Compute overlap of alignments on the SVs
resultsGroupped['coverage'] = [resultsGroupped.loc[i,'svBED'].overlapLen(resultsGroupped.loc[i,'alnBED'], percent = True) for i in resultsGroupped.index]

# Keep only SVs covered on more than 50%
resultsGroupped = resultsGroupped.loc[resultsGroupped['coverage'] >= 50,:]

# Output IDs of Ty related SVs
resultsGroupped.loc[:,'qseqid'].to_csv(outputPath, sep = '\t', index = False, header = False)

# Cleaning ==================================================================
cmd = ['rm', '-f', 
       'TyDB.ndb', 
       'TyDB.nhr', 
       'TyDB.nin', 
       'TyDB.njs', 
       'TyDB.not', 
       'TyDB.nsq', 
       'TyDB.ntf', 
       'TyDB.nto',
       'TySearch.blastn']
subprocess.call(cmd)
