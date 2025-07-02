#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script retrieve contigs longer than X kb. 

It takes as input :
	-f --fasta: a fasta file (multi fasta)
	-m --minLength: minimum length of a contig (kb)
	-o --output: output prefix


blast+ is used in this script. If blast+ is not in the path, the path can
be added to the variable blast. 
'''
# ---------------------------------------------------------------------------
import sys
import argparse
from Tools import *
# ---------------------------------------------------------------------------
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(
	description = """This script retrieve contigs longer than X kb. """)
parser.add_argument("-f", "--fasta", help="Fasta file (multi fasta)", required=True)
parser.add_argument("-m", "--minLength", help="minimum length of a contig (kb)", type=float, default=1)
parser.add_argument("-o", "--output", help="Prefix of the output file", default = "")

# Read arguments from the command line
args = parser.parse_args()

fastaPath=args.fasta
minLength=args.minLength
outputPath=args.output

# ===============
# Get Input files
# ===============
# Read fasta
fasta = Fasta(fastaPath)

# Check each sequence in Fasta
toKeep = Fasta()
for seq in fasta:
	if len(seq) >= minLength * 1000:
		toKeep += Fasta([seq])

# write output file
if outputPath != "":
	toKeep.toFile(f'{outputPath}.min{minLength}kb.fasta')
else:
		sys.stdout.write(toKeep.__str__()+'\n')
