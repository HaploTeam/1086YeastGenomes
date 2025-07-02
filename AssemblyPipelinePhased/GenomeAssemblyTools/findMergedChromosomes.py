#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script check if there are contigs corresponding to merged chromosomes
of a reference. It looks if a single contigs aligns on more than 90% of 2 
reference chromosomes. 
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentToMatch: Percentage a contig has to cover on each chromosome to be considered as a merge
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be:
Assembly	NbMergedContigs
Draft1		0
Draft2		0
Draft3		1
'''
# ---------------------------------------------------------------------------
import os
import argparse
import time
from datetime import datetime
from random import randint
import csv
from Tools import *
# ---------------------------------------------------------------------------
# Definitions
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
parser = argparse.ArgumentParser(description = """This script check if there are contigs corresponding to merged chromosomes
of a reference. It looks if a single contigs aligns on more than 90% of 2 
reference chromosomes. 
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions
""")
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Name of the output file", required=True)
parser.add_argument("-p", "--percentToMatch", help="Percentage a contig has to cover on each chromosome to be considered as a merge", type=int, default=90)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=20)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
percent=args.percentToMatch
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

# Write header in output file
out = open(outputPath, 'w')
out.write("Assembly\tNbMergedContigs\n")

# ========================================
# Get reference chromosome name and length
refFasta = Fasta(refPath)
refName=refPath.split("/")[-1]

# ==================================================
# Iterate over each draft assembly to get the number 
# of contigs that cover X% of the reference
# ==================================================
for d in range(nbDraft):

	draftName=draftPaths[d].split("/")[-1]
	print("\n\t--- Running for draft assembly "+draftName+" ---\n")
	startTime = time.time()

	# Get draft contig name
	draftFasta = Fasta(draftPaths[d])

	# Align draft to ref
	prefix = "Alignment_" + refName + "_" + draftName + "_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))
	getShowCoords(refPath, draftPaths[d], prefix, mummer, threads)
	endTime = time.time()
	print("1/2 Alignment done: ran in "+str(round(endTime-startTime))+"s")

	startTime = time.time()
	# Get BED list of reference genome
	refBED = getBED(refPath)

	# Create a matrix which will contain BEDcoordinates objects of covered regions
	# First dimension: draft contigs
	# Second dimension: reference chromosomes
	x = len(draftFasta)
	y = len(refFasta)
	analysisMatrix = []
	for i in range(x):
		analysisMatrix.append([])
		for j in range(y):
			analysisMatrix[i].append([])

	# For each line in the alignment coords file, add BEDcoordinates of each alignment
	coordsFile = open(prefix+".coords")
	coords = csv.reader(coordsFile, delimiter="\t")
	for row in coords:
		draftContigIndex=draftFasta.getIndexFromID(row[8])
		refChrIndex=refFasta.getIndexFromID(row[7])
		# get BEDcoordinates to add
		chr=row[7]
		startPos=int(row[0])
		endPos=int(row[1])
		coord2add = BEDcoordinates(id = chr, start = startPos, end = endPos)
		# Add BEDcoordinates
		analysisMatrix[draftContigIndex][refChrIndex] += [coord2add]
	coordsFile.close()
	os.remove(prefix+".coords")

	# Merge all BEDcoordinates in one BED per point in the matrix
	for i in range(len(draftFasta)):
		for j in range(len(refFasta)):
			analysisMatrix[i][j] = BED(analysisMatrix[i][j])

	# For each point of the matrix, get the percentage of reference chromosome covered
	for i in range(len(draftFasta)):
		for j in range(len(refFasta)):
			# Get ref chr BED
			refChrBED = refBED.getID(refFasta.getID()[j])
			# Get proportion of chromosome covered
			analysisMatrix[i][j] = refChrBED.overlapLen(analysisMatrix[i][j], percent = True)

	# For each draft contig, check if it is not 2 merged reference chromosomes
	nbContigsMerged = 0
	for i in range(len(draftFasta)):
		# If coverage is superior to 90% for more than one ref chromosome: it is a merged contig
		supToX = len([x for x in analysisMatrix[i] if x >= percent])
		if supToX > 1:
			nbContigsMerged += 1

	endTime = time.time()
	print("2/2 Coordinates analsis done: ran in "+str(round(endTime-startTime))+"s\n")

	out.write(draftName + "\t" + str(nbContigsMerged) + "\n")
	print("Assembly\tNbMergedContigs")
	print(draftName + "\t" + str(nbContigsMerged))

out.close()




