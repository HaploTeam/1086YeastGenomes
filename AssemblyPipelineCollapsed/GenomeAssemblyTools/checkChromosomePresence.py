#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script check if all chromosomes of a reference are covered by
a draft assembly. 
It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta, several allowed)
	-o --output: name of output tsv file
	-p --percentCovered: Percentage of each chromosome that have to be covered
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer
	-i --detailledInfo: Display the coverage of each chromosome of the reference

Output will be:
Assembly	NbChrCovered
Draft1		16
Draft2		15
Draft3		16
'''
# ---------------------------------------------------------------------------
import os
import sys
import argparse
import time
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
parser = argparse.ArgumentParser(
	description = """This script check if all chromosomes of a reference are covered by a draft assembly. 
	It uses the nucmer (--maxmatch) and show-coords MUMmer4's functions. 
	The percentage of each chromosome to be covered can be specified for all chromosomes, 
	or each one specifically. """)
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", nargs='+', required=True)
parser.add_argument("-o", "--output", help="Output file", default = "")
parser.add_argument("-p", "--percentCovered", help="Percentage of each chromosome that have to be covered", type=int, default=80)
parser.add_argument("-pid", "--percentCoveredForID", help="Percentage that have to be covered for a specific chromosome. Syntax is chromosome1=75 for min coverage of 75p on chromosome1. Several allowed after -pid. ", nargs='+', default = "", type=str)
parser.add_argument("-i", "--detailledInfo", help="Display the coverage of each chromosome of the reference", action='store_true')
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=20)

# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPaths=args.draft
nbDraft=len(draftPaths)
outputPath=args.output
percent=int(args.percentCovered)
percentID=args.percentCoveredForID
nbIDspecific = len(percentID)
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

# Write header in output file
if args.detailledInfo:
	if outputPath != "":
		out = open(outputPath, 'w')
		out.write("Assembly\tChromosome\tCoverage\n")
	else:
		sys.stdout.write("Assembly\tChromosome\tCoverage\n")
else:
	if outputPath != "":
		out = open(outputPath, 'w')
		out.write("Assembly\tNbChrPresent\n")
	else:
		sys.stdout.write("Assembly\tNbChrPresent\n")

# ========================================
# Get reference chromosome name and length
# Read reference fasta
refFasta = Fasta(refPath)

# make list of percentage to cover for each chromosome
percentList = [percent] * len(refFasta) # Default percentage for every Chr
if nbIDspecific != 0: # If specific percentage for a chromosome
	for i in range(nbIDspecific):
		specChr = percentID[i].split("=")[0]
		specPerc = int(percentID[i].split("=")[1])
		if specChr not in refFasta.getID():
			raise ValueError("Chromosome ID specified for specific coverage is not found in the reference genome. ")
		else:
			percentList[refFasta.getIndexFromID(specChr)] = specPerc

refName=refPath.split("/")[-1]

# ==================================================
# Iterate over each draft assembly to get the number 
# of chromosomes covered by the draft assembly
# ==================================================
for d in range(nbDraft):

	draftName=draftPaths[d].split("/")[-1]
	#print("\n\t--- Running for draft assembly "+draftName+" ---\n")
	start = time.time()

	# Align draft to ref
	prefix = "Alignment_" + refName + "_" + draftName + "_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))
	getShowCoords(refPath, draftPaths[d], prefix, mummer, threads)
	end = time.time()
	#print("1/2 Alignment done: ran in "+str(round(end-start))+"s")

	start = time.time()
	# Get BED list of reference genome
	refBED = getBED(refPath)

	# ======================================
	# Get BED of alignments for each chromosome of the ref
	alignmentBEDs = [[] for i in refFasta]
	# Alignments BED will contain per Chromosome the BED of alignments of this chromosome on the draft genome
	with open(prefix+".coords", 'r') as coords:
		for line in coords:
			refStart = int(line.split("\t")[0])
			refEnd = int(line.split("\t")[1]) + 1
			Chr = line.split("\t")[7]
			ChrIndex = refFasta.getIndexFromID(Chr)

			alignmentBEDs[ChrIndex] += [BEDcoordinates(id = Chr, start = refStart, end = refEnd)]
	
	os.remove(prefix+".coords")

	# Convert all to BED
	alignmentBEDs = [BED(x) for x in alignmentBEDs]

	# Get coverage per chromosome
	refCoverage = []
	for i in range(len(refFasta)):
		refCoverage += [refBED.getID(refFasta.getID()[i]).overlapLen(alignmentBEDs[i], percent = True)]

	# Output detailled for each chromosome
	if args.detailledInfo:
		for i in range(len(refFasta)):
			if outputPath != "":
				out.write(draftName + "\t" + refFasta.getID()[i] + "\t" + str(refCoverage[i]) + "\n")
			else: 
				sys.stdout.write(draftName + "\t" + refFasta.getID()[i] + "\t" + str(refCoverage[i]) + "\n")

	# Output total number of chromosome presents
	else:
		# Get number of chromosome covered at X%
		nbChrPresent = 0
		for i in range(len(refFasta)):
			if refCoverage[i] >= percentList[i]:
				nbChrPresent += 1

		if outputPath != "":
			out.write(draftName + "\t" + str(nbChrPresent) + "\n")
		else: 
			sys.stdout.write(draftName + "\t" + str(nbChrPresent) + "\n")

if outputPath != "":
	out.close()




