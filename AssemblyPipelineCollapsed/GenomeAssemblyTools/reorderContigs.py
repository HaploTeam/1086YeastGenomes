#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script reorder a fasta file according to a reference fasta file

It takes as input :
	-r --ref: the reference genome assembly (multi fasta)
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-o --output: output file
	-t --threads: Number of threads to use for Mummer, optional (default = 1)
	-m --mummerPath: Path to mummer function, if not in path

MUMmer4 is used in this script. If MUMmer4 is not in the path, the path can
be added to the variable mummer. 
'''
# ---------------------------------------------------------------------------
import csv
import os
import argparse
from datetime import datetime
from random import randint
from Tools import *
# ---------------------------------------------------------------------------
# Definitions

# def reverse_complement(dna):
# 	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N', 'S':'S', 'W':'W', 'Y':'R', 'R':'Y', 'M':'K', 'K':'M', 'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'n':'n', 's':'s', 'w':'w', 'y':'r', 'r':'y', 'm':'k', 'k':'m'}
# 	return ''.join([complement[base] for base in dna[::-1]])

def rank_simple(vector):
	return sorted(range(len(vector)), key=vector.__getitem__)

# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = 
'''
This script reorder a fasta file according to a reference fasta file. 
MUMmer4 is used in this script. If MUMmer4 is not in the path, the path can
be added to the variable mummer. 
'''
)
parser.add_argument("-r", "--ref", help="reference genome assembly (multi fasta)", required=True)
parser.add_argument("-d", "--draft", help="draft genome assembly (multi fasta)", required=True)
parser.add_argument("-o", "--output", help="Name of the outputed fasta", required=True)
parser.add_argument("-t", "--threads", help="Number of threads to use for Mummer (default = 1)", type=int, default=1)
parser.add_argument("-m", "--mummerPath", help="Path to mummer function, if not in path", type=str, default="")

# Read arguments from the command line
args = parser.parse_args()

refPath=args.ref
draftPath=args.draft
outputPath=args.output
nbThreads=args.threads
mummer=args.mummerPath
if mummer != "" and not mummer.endswith("/") :
	mummer += "/"

print("\n\t--- REORDERING CONTIGS ---\n")
print("Arguments detected:")
print("\t--ref:\t"+refPath)
print("\t--draft:\t"+draftPath)
print("\t--output:\t"+outputPath)
if nbThreads != 1:
	print("\t--threads:\t"+str(nbThreads))
print('')

# ===============
# Get Input files
# ===============

# read reference Fasta
print("1/7\tGet reference chromosomes names")
refFasta = Fasta(refPath)

# Get draft contigs names and sequences
print("2/7\tGet draft contigs names and sequences")
draftFasta = Fasta(draftPath)

# =========================================
# Using MUMmer4 to align draft to reference
# =========================================

# Align draft to reference using nucmer
print("3/7\tAlign draft to reference using nucmer")
prefix = "Alignment_" + datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000))

nucmerShellCommand=mummer+"nucmer --prefix "+prefix+" -t "+str(nbThreads)+" "+refPath+" "+draftPath
os.system(nucmerShellCommand)

# Get alignment coordinates
print("4/7\tGet alignment coordinates")
showcoordsShellCommand=mummer+"show-coords -THq "+prefix+".delta > "+prefix+".coords"
os.system(showcoordsShellCommand)

# ==================
# Alignment analysis
# ==================

# Read alignment coordinates
print('5/7\tRead alignment coordinates')
# For each draft contig, get the total alignment length on each ref chromosome, the mean position on reference chromosome and the number of forward and reverse alignments
alignmentLength=[] 		# 0 # total length of alignment of the draft contig on each ref chromosome
alignmentPos=[] 		# 1 # Middle of the longest alignment. tuple (middle position, alignment length)
alignmentNb=[] 			# 2 # Nb of alignments on each ref Chr
alignmentNbForward=[] 	# 3 # Length of forward alignments on each ref Chr
alignmentNbReverse=[] 	# 4 # Length of reverse alignments on each ref Chr
# 5 statistics are searched

# Create a 3D matrix containing all statistics
# First level: draft contig
# Second level: type of statistics (alignmentLength, alignmentPos, alignmentNb, alignmentNbReverse, alignmentNbForward)
# Third level: ref chromosome
x = len(draftFasta)
y = 5
z = len(refFasta)
analysisMatrix = []
for i in range(x):
	analysisMatrix.append([])
	for j in range(y):
		analysisMatrix[i].append([])
		for k in range(z):
			if j == 1:
				# Add tuple for middle position and length of the longuest alignment
				analysisMatrix[i][j].append((0,0))
			else:
				analysisMatrix[i][j].append(0)


with open(prefix+".coords") as coordsFile:
	coords = csv.reader(coordsFile, delimiter="\t")
	for row in coords:
		draftContigIndex=draftFasta.getIndexFromID(row[8])
		refChrIndex=refFasta.getIndexFromID(row[7])

		# Add alignment length
		analysisMatrix[draftContigIndex][0][refChrIndex] += int(row[4])
		# Add alignment middle position if it is the longuest alignment
		if analysisMatrix[draftContigIndex][1][refChrIndex][1] < int(row[4]):
			analysisMatrix[draftContigIndex][1][refChrIndex] = ( int( (int(row[0]) + int(row[1])) / 2 ), int(row[4]))
		# Add number of alignments
		analysisMatrix[draftContigIndex][2][refChrIndex] += 1
		# Add length of reverse or foward alignments
		if int(row[2]) < int(row[3]):
			analysisMatrix[draftContigIndex][3][refChrIndex] += int(row[4])
		else:
			analysisMatrix[draftContigIndex][4][refChrIndex] += int(row[4])



# For each draft contig, get the corresponding ref Chr, orientation, and position
print("6/7\tFind corresponding chromosomes")
correspondingChrs=[""]*len(draftFasta)
orientations=[""]*len(draftFasta)
positions=[""]*len(draftFasta)

for draftContigIndex in range(len(draftFasta)):
	# Find the corresponding chromosome on the reference
	maxLength = max(analysisMatrix[draftContigIndex][0])
	refChrIndex = analysisMatrix[draftContigIndex][0].index(maxLength)
	correspondingChr = refFasta.getID()[refChrIndex]
	correspondingChrs[draftContigIndex] = correspondingChr

	# Find the orientation of the contig
	forwardLength = analysisMatrix[draftContigIndex][3][refChrIndex]
	reverseLength = analysisMatrix[draftContigIndex][4][refChrIndex]
	if forwardLength > reverseLength:
		orientation = "F"
	else:
		orientation = "R"
	orientations[draftContigIndex] = orientation

	# get the middle position of the longuest alignment of the contig on the ref chr
	position=analysisMatrix[draftContigIndex][1][refChrIndex][0]
	positions[draftContigIndex] = position

print("\nCorresponding chromosomes")
for i in range(len(draftFasta)):
	print(f'{draftFasta.getID()[i]}\t{correspondingChrs[i]}')

# get the final draft contigs order
order=[0]*len(draftFasta)

index=0
for chr in refFasta.getID():
	if correspondingChrs.count(chr) == 0:
		pass
	elif correspondingChrs.count(chr) == 1:
		draftContigIndex = correspondingChrs.index(chr)
		order[draftContigIndex] = index
		index += 1
	else:
		# Get indices of all contigs aligning on the same chromosome
		draftContigIndices = [i for i, x in enumerate(correspondingChrs) if x == chr]
		# Get contigs id
		contigs = [draftFasta.getID()[i] for i in draftContigIndices]
		# Get relative position of contigs on the reference chromosome
		posOnChr = [positions[i] for i in draftContigIndices]
		# Get ranks of contigs
		innerOrder = rank_simple(posOnChr)
		# Sort contigs by position order
		contigsOrdered = [contigs[i] for i in innerOrder]
		# Add final order
		for c in contigsOrdered:
			draftContigIndex = draftFasta.getIndexFromID(c)
			order[draftContigIndex] = index
			index += 1


# Get the list of ordered index (reverse of order)
orderToPrint=[]
for i in range(len(order)):
	orderToPrint += [order.index(i)]

# =======================
# Write output fasta file
# =======================

sequences = []
for index in orderToPrint:
	seq = draftFasta.sequences[index]
	if orientations[index] == "R":
		seq.reverseComplement()
	sequences += [seq]

newDraftFasta = Fasta(sequences)

print("\n7/7\tWrite to output file: "+outputPath)
newDraftFasta.toFile(outputPath)

# =====================
# clean temporary files
# =====================

os.system("rm -f "+prefix+".coords")
os.system("rm -f "+prefix+".delta")

print("\n\t--- REORDERING RAN SUCCESFULLY ---\n")

