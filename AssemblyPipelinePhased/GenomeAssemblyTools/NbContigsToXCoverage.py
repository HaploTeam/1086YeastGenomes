#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script check how many contigs are needed to cover X% of the reference
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentCovered: Percentage of the genome which has to be covered
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be:
Assembly	NbContigsToX
Draft1		16
Draft2		17
Draft3		16
'''
# ---------------------------------------------------------------------------
import os
import argparse
import time
import sys
from datetime import datetime
from random import randint
from Tools import *
# ---------------------------------------------------------------------------
# Definitions

def rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)

def decreasing_rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)[::-1]

def getShowCoords(refPath, draftPath, prefix, mummerPath, threads):
	# Run nucmer
	nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --maxmatch --prefix "+prefix+" "+refPath+" "+draftPath
	os.system(nucmerShellCommand)
	# Run show coords
	showcoordsShellCommand=mummerPath+"show-coords -TH "+prefix+".delta > "+prefix+".coords"
	os.system(showcoordsShellCommand)
	# Remove nucmer delta file
	os.remove(prefix+".delta")

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = 
'''
This script check how many contigs are needed to cover X% of the reference
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions. 
'''
)
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Output file", default = "")
parser.add_argument("-p", "--percentCovered", help="Percentage of the genome which has to be covered", type=int, default=95)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=20)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
percent=args.percentCovered
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

# Write header in output file
if outputPath != "":
	out = open(outputPath, 'w')
	out.write("Assembly\tNbContigsTo" + str(percent) + "\n")
else:
	print("Assembly\tNbContigsTo" + str(percent))

# ========================================
# Get ref file name
refName=refPath.split("/")[-1]

# ==================================================
# Iterate over each draft assembly to get the number 
# of contigs that cover X% of the reference
# ==================================================
for d in range(nbDraft):

	draftName=draftPaths[d].split("/")[-1]
	#print("\n\t--- Running for draft assembly "+draftName+" ---\n")
	start = time.time()

	# Get draft contig name
	draftChr=[]
	with open(draftPaths[d], 'r') as draft:
		for line in draft:
			if line.startswith(">"):
				draftChr += [line.split(">")[1].split()[0]]

	# Align draft to ref
	prefix = "Alignment_" + refName + "_" + draftName + "_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))
	getShowCoords(refPath, draftPaths[d], prefix, mummer, threads)
	end = time.time()
	#print("1/2 Alignment done: ran in "+str(round(end-start))+"s")


	start = time.time()
	# Get BED list of reference genome
	refBED = getBED(refPath)

	# ======================================
	# Get reference coverage for each contig
	alignmentBEDs = [[] for i in range(len(draftChr))]
	# Alignments BED will contain per contig the BED of alignments of this contig on the reference genome
	with open(prefix+".coords", 'r') as coords:
		for line in coords:
			refStart = int(line.split("\t")[0])
			refEnd = int(line.split("\t")[1]) + 1
			refChr = line.split("\t")[7]
			contig = line.strip().split("\t")[8]
			contigIndex = draftChr.index(contig)

			alignmentBEDs[contigIndex] += [BEDcoordinates(id = refChr, start = refStart, end = refEnd)]

	# Convert all to BED
	alignmentBEDs = [BED(x) for x in alignmentBEDs]

	# Get coverage per contig
	draftCoverage = []
	for i in range(len(draftChr)):
		draftCoverage += [refBED.overlapLen(alignmentBEDs[i], percent = True)]

	# Order contigs by decreasing coverage of reference assembly
	order = decreasing_rank_simple(draftCoverage)

	# =========================================================
	# Get reference coverage by adding sequentially each contig

	contigs = []
	coords = [] # coords will contain all the BEDcoordinates corresponding to all the contigs
	increasingDraftCoverage = []
	for i in order:
		contigs += [i]
		coords += [alignmentBEDs[contigs[-1]]] # Add BEDcoords corresponding to last contig added
		c = refBED.overlapLen(BED(coords), percent = True)
		increasingDraftCoverage += [c]

	nbContigs = 0
	for c in increasingDraftCoverage:
		nbContigs += 1
		if c >= percent:
			break
	if c >= percent:
		n = nbContigs
	else:
		n = "NotReached"
	os.remove(prefix+".coords")

	if outputPath != "":
		out.write(f'{draftName}\t{n}\n')
	else:
		sys.stdout.write(f'{draftName}\t{n}\n')

if outputPath != "":
	out.close()




