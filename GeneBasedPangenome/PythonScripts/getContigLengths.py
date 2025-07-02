#!/usr/bin/env python
# -*- coding: utf8 -*-
# This script output length of sequences in a fasta file
import argparse
import sys
sys.path.insert(1, '/home/vloegler/Pangenome_Sace1000ONT/02-scripts/GenomeAssemblyTools')
from Tools import *
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = """
This script output length of sequences in a fasta file. 
""")
parser.add_argument("-f", "--fasta", help="fasta file", required=True)

# Read arguments from the command line
args = parser.parse_args()

fastaPath = args.fasta
fasta = Fasta(fastaPath)

for s in fasta.sequences:
    sys.stdout.write(f'{s.id}\t{len(s)}\n')


