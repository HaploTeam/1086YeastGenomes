#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2023/10/04
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script changes the ID of the sequences in a fasta file. It requires a
tab-separated file with col1: old ID and col2: new ID. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import argparse
import os
import logging
import re
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
parser.add_argument("-id", "--id", help="Tab-separated file with old IDs (col1) and new IDs (col2). No header. ", required=True)
parser.add_argument("-o", "--output", help="Output file (fasta format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

fastaPath=os.path.abspath(args.fasta)
idPath=os.path.abspath(args.id)
outputPath=os.path.abspath(args.output)

#fastaPath = 'SV.1087Samples.InsertionSequences.fasta'
#idPath='IDs.txt'

# Read ID table
IDs = pd.read_csv(idPath, sep = '\t', names = ['OldID', 'NewID'])
IDs.index = IDs['OldID']

# Read Fasta
fasta = Fasta(fastaPath)

# Create df with fasta sequences
fastaDF = pd.DataFrame({'Obj': fasta.sequences})

# Extract ID, header and Sequence
fastaDF['ID'] = [s.id for s in fastaDF['Obj']]
fastaDF.index = fastaDF['ID']
fastaDF['seq'] = [s.seq for s in fastaDF['Obj']]
fastaDF['description'] = [s.description for s in fastaDF['Obj']]

# Remove ID equivalences that are not in the fasta file
IDs = IDs.loc[IDs['OldID'].isin(fastaDF['ID'])]

# Add New ID column and change IDs that are specified in the tab file
fastaDF['NewID'] = fastaDF['ID']
fastaDF.loc[IDs['OldID'],'NewID'] = IDs['NewID']
logging.info(f'{IDs.shape[0]} sequence IDs changed. ')

# Replace ID in description column
fastaDF['NewDescription'] = [re.sub(fastaDF.loc[i, 'ID'], fastaDF.loc[i, 'NewID'], fastaDF.loc[i, 'description']) for i in fastaDF.index]

# Create new Sequence objects with the new description
fastaDF['NewObj'] = [Sequence(fastaDF.loc[i, 'NewDescription'], fastaDF.loc[i, 'seq']) for i in fastaDF.index]

# Regenerate fasta object with changes IDs and write output file
newFasta = Fasta(fastaDF['NewObj'].tolist())
newFasta.toFile(outputPath)


