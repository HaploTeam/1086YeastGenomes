#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/28
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script rename contigs according to a reference genome, if at least X kb 
or X% of a contig aligns on a reference chromosome. If a contig aligns on 
several reference chromosomes, all chromosomes name will be indicated, by 
order of appearence in the contig. 
An additionnal argument add the alignment sizes in the contig name. 
This script uses nucmer and show-coords MUMmer4's functions. 

It's better to use a reference with transposable elements masked, for example
with RepeatMasker:
	RepeatMasker reference.fasta -lib TransposableElements.fasta -s -e ncbi
and with telomeres masked.

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: draft genome assemblies (multi fasta)
	-o --output: name of output fasta file
	-k --alignmentSizeKb: Minimum alignment size to have name correspondance (in kb)
	-p --alignmentSizePercent: Minimum alignment size to have name correspondance (in % of contig size)
	-P --prefix: Prefix to add before the chromosome name
	-i --alignmentSizeInfo: Add alignment size information in contig names
	-m --mummerPath: Path to mummer function, if not in path
	-t --threads: Number of threads for nucmer

Output will be fasta with renamed contigs. 
'''
# ---------------------------------------------------------------------------
import os
import argparse
import csv
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
	#nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --maxmatch --prefix "+prefix+" "+refPath+" "+draftPath
	nucmerShellCommand=mummerPath+"nucmer -t "+str(threads)+" --mum --prefix "+prefix+" "+refPath+" "+draftPath
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
This script rename contigs according to a reference genome, if at least X kb 
or X% of a contig aligns on a reference chromosome. If a contig aligns on 
several reference chromosomes, all chromosomes name will be indicated, by 
order of appearence in the contig. 
An additionnal argument add the alignment sizes in the contig name. 
This script uses nucmer and show-coords MUMmer4's functions. 

It's better to use a reference with transposable elements masked, for example
with RepeatMasker:
	[RepeatMasker reference.fasta -lib TransposableElements.fasta -s -e ncbi]
and with telomeres masked.

Output will be fasta with renamed contigs. 
'''
)
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assemblies (multi fasta)", required=True)
parser.add_argument("-o", "--output", help="Name of the outputed Fasta", required=True)
parser.add_argument("-k", "--alignmentSizeKb", help="Minimum alignment size to have name correspondance (in kb, default 50)", type=int, default=50)
parser.add_argument("-p", "--alignmentSizePercent", help="Minimum alignment size to have name correspondance (in percent, default 20)", type=int, default=20)
parser.add_argument("-P", "--prefix", help="Prefix to add before the chromosome name", type=str, default="")
parser.add_argument("-i", "--alignmentSizeInfo", help="Add alignment size information in contig names", action='store_true')
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")
parser.add_argument("-t", "--threads", help="Number of threads for nucmer", type=int, default=20)


# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPath=args.draft
outputPath=args.output
alignmentSizeKb=args.alignmentSizeKb
alignmentSizePercent=args.alignmentSizePercent
ChrPrefix=args.prefix
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"
threads=args.threads

print("\n\t--- RENAMING CONTIGS ---\n")
print("Arguments detected:")
print("\t--ref:\t"+refPath)
print("\t--draft:\t"+draftPath)
print("\t--output:\t"+outputPath)
print("\t--alignmentSizeKb:\t"+str(alignmentSizeKb)+"kb")
print("\t--alignmentSizePercent:\t"+str(alignmentSizePercent)+"%")
if ChrPrefix != "":
	print("\t--prefix:\t"+ChrPrefix)
if args.alignmentSizeInfo:
	print("\t--alignmentSizeInfo")
if mummer != "":
	print("\t--mummerPath:\t"+mummer)
if threads != 1:
	print("\t--threads:\t"+str(threads))
print('')

# =============================
# Read reference fasta
refFasta = Fasta(refPath)

# =================================
# Read draft fasta
draftFasta = Fasta(draftPath)


refName=refPath.split("/")[-1]
draftName=draftPath.split("/")[-1]
prefix = "Alignment_" + refName + "_" + draftName + "_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))

# ===============================================
# Align and get alignment coordinates with nucmer
getShowCoords(refPath, draftPath, prefix, mummer, threads)

# ==========================
# Read alignment coordinates
# Cummulative size of alignments will be stored in a matrix with 
# first dimension: Draft contig
# second dimension: Reference chromosome
# analysisMatrixRef contains BED based on ref chr
# analysisMatrixDraft contains BED based on draft chr
x = len(draftFasta)
y = len(refFasta)
analysisMatrixRef = []
analysisMatrixDraft = []
for i in range(x):
	analysisMatrixRef.append([])
	analysisMatrixDraft.append([])
	for j in range(y):
		analysisMatrixRef[i].append([])
		analysisMatrixDraft[i].append([])

# Fill matrix
with open(prefix+".coords") as coordsFile:
	coords = csv.reader(coordsFile, delimiter="\t")
	for row in coords:
		draftContigIndex=draftFasta.getIndexFromID(row[8])
		refChrIndex=refFasta.getIndexFromID(row[7])

		# Add alignment length to analysisMatrixRef
		chr=row[7]
		startPos=int(row[0])
		endPos=int(row[1]) + 1
		if startPos < endPos:
			coord2add = BEDcoordinates(id = chr, start = startPos, end = endPos)
		else:
			coord2add = BEDcoordinates(id = chr, start = endPos, end = startPos)
		analysisMatrixRef[draftContigIndex][refChrIndex] += [coord2add]

		# Add alignment length to analysisMatrixDraft
		chr=row[8]
		startPos=int(row[2])
		endPos=int(row[3]) + 1
		if startPos < endPos:
			coord2add = BEDcoordinates(id = chr, start = startPos, end = endPos)
		else:
			coord2add = BEDcoordinates(id = chr, start = endPos, end = startPos)
		analysisMatrixDraft[draftContigIndex][refChrIndex] += [coord2add]

os.remove(prefix+".coords")

# Convert all point of the matrix to BED object and get length
for i in range(len(draftFasta)):
	for j in range(len(refFasta)):
		analysisMatrixDraft[i][j] = BED(analysisMatrixDraft[i][j])
		analysisMatrixRef[i][j] = BED(analysisMatrixRef[i][j])

# Second matrix just containing length
x = len(draftFasta)
y = len(refFasta)
analysisMatrixLen = []
for i in range(x):
	analysisMatrixLen.append([])
	for j in range(y):
		analysisMatrixLen[i].append(analysisMatrixRef[i][j].len)

# Find correspondance between contigs and chromosomes
# A correspondance is an total alignment size greater than 30kb, 
# 20% the contig size or 20% the chromosome size
# Results will be stored in a matrix correspondanceMatrixRef[refIndex][draftIndex]
# that will contain either 0 (no correspondance between contig and chr)
# or a number which is the part of Chromosome corresponding
# 1 = first part or whole chromosome
# 2 = second part
# ...
#
# A second matrix correspondanceMatrixDraft will contain the order in which the 
# chromosome fragments are present in the contig

x = len(draftFasta)
y = len(refFasta)
correspondanceMatrixRef = []
correspondanceMatrixDraft = []
for i in range(x):
	correspondanceMatrixRef.append([])
	correspondanceMatrixDraft.append([])
	for j in range(y):
		correspondanceMatrixRef[i].append(0)
		correspondanceMatrixDraft[i].append(0)

# Find fragment of chromosomes
#numberFragments = [0]*len(refChr)
for j in range(len(refFasta)):
	contigAssociated = [0]*len(draftFasta)
	for i in range(len(draftFasta)):
		sizeThreshold = min([alignmentSizeKb*1000, alignmentSizePercent*0.01*refFasta.getLengths()[j], alignmentSizePercent*0.01*draftFasta.getLengths()[i]])
		if analysisMatrixLen[i][j] >= sizeThreshold:
			contigAssociated[i] = analysisMatrixRef[i][j].getCenter()

	# Order chromosome fragments
	contigAssociatedOrder = [0]*len(draftFasta)
	centers = [x for x in contigAssociated if x > 0]
	for i in range(len(draftFasta)):
		if contigAssociated[i] != 0:
			correspondanceMatrixRef[i][j] = 1 + sorted(centers).index(contigAssociated[i])



for i in range(len(draftFasta)):
	chromoAssociated = [0]*len(refFasta)
	for j in range(len(refFasta)):
		sizeThreshold = min([alignmentSizeKb*1000, alignmentSizePercent*0.01*refFasta.getLengths()[j], alignmentSizePercent*0.01*draftFasta.getLengths()[i]])
		if analysisMatrixLen[i][j] >= sizeThreshold:
			chromoAssociated[j] = analysisMatrixDraft[i][j].getCenter()

	# Order chromosome in contigs
	chromoAssociatedOrder = [0]*len(refFasta)
	centers = [x for x in chromoAssociated if x > 0]
	for j in range(len(refFasta)):
		if chromoAssociated[j] != 0:
			correspondanceMatrixDraft[i][j] = 1 + sorted(centers).index(chromoAssociated[j])

# Get name of each contig
draftNewChr = []
unplacedNb=0
for i in range(len(draftFasta)):
	newName = ">"
	if ChrPrefix != "":
		newName += ChrPrefix+"_"

	additionnalInfo = " AlignmentSizeOnRef:"
	nbChrInContig = max(correspondanceMatrixDraft[i]) # Number of chr in the contig
	if nbChrInContig > 0:
		for n in range(nbChrInContig):
			if n != 0:
				newName += "_"
				additionnalInfo += ","
			# Add Chr Name
			j = correspondanceMatrixDraft[i].index(n+1)
			newName += refFasta.getID()[j]
			# Add Chr fragment, if chromosome if fragmented
			if max([x[j] for x in correspondanceMatrixRef]) > 1:
				newName += "."+str(correspondanceMatrixRef[i][j])
			# Add alignment size
			additionnalInfo += refFasta.getID()[j] + "=" + str(analysisMatrixLen[i][j])
	else:
		unplacedNb += 1
		newName += "Unplaced." + str(unplacedNb)
		additionnalInfo = ""

	if args.alignmentSizeInfo:
		draftNewChr += [newName+" length="+str(draftFasta.getLengths()[i])+additionnalInfo]
	else:
		draftNewChr += [newName+" length="+str(draftFasta.getLengths()[i])]

# Change ID in a new draft fasta object
newDraftFasta = Fasta()
index = 0
for seq in draftFasta:
	# Change ID
	seq.id = draftNewChr[index].split()[0]
	# Change sequence description
	seq.description = f'{draftNewChr[index]}\n'
	# Add to new fasta object
	newDraftFasta += Fasta([seq])
	index += 1


# Write new fasta object to file
newDraftFasta.toFile(outputPath)

print(f"New fasta written to: {outputPath}\n")

