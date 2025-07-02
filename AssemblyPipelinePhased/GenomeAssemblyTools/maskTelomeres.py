#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/02/24
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script masks telomeres (specified by a size) of a given assembly. 
Nucleotides are replaced by N. 

It takes as input :
	-f --fasta: fasta file
	-o --output: output file
	-k --kbToHide: Size of regions to hide at chromosome extremities(kb, default = 20)
'''
# ---------------------------------------------------------------------------
import sys
import argparse
from Tools import *
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = 
	'''
	This script masks telomeres (specified by a size) of a given assembly. 
	Nucleotides are replaced by N. 
	'''
	)
parser.add_argument("-f", "--fasta", help="fasta file", required=True)
parser.add_argument("-o", "--output", help="output file", default = "")
parser.add_argument("-k", "--kbToHide", help="Size of regions to hide at chromosome extremities(kb, default = 20)", type=float, default=20)

# Read arguments from the command line
args = parser.parse_args()

fastaPath = args.fasta
outputPath=args.output
bpToHide=int(args.kbToHide*1000)

# ===============
# Get Input files
# ===============

# Read fasta
fasta = Fasta(fastaPath)

for s in fasta:
	s.seq = "N"*bpToHide + s.seq[bpToHide:len(s)-bpToHide] + "N"*bpToHide

if outputPath != "":
	fasta.toFile(outputPath)
else:
	sys.stdout.write(fasta.__str__()+'\n')
