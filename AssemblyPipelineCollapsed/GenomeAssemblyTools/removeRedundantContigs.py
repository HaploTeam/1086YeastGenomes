#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script remove redundant contigs. If one contig is fully covered 
(95% threshold) by other contigs of the draft assembly, the contig is removed. 

It takes as input :
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-o --output: output prefix
	-b --blastPath: Path to blast+, if not in path


blast+ is used in this script. If blast+ is not in the path, the path can
be added to the variable blast. 
'''
# ---------------------------------------------------------------------------
import csv
import os
import argparse
from datetime import datetime
from random import randint
import time
from Tools import *
# ---------------------------------------------------------------------------
# Definitions
def blastn(fasta1, fasta2, blastPath, out):
	blastCommand = blastPath + "blastn -query " + fasta1 + " -subject " + fasta2 + " -outfmt 6 -out " + out
	print(blastCommand)
	os.system(blastCommand)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = 
'''
This script remove redundant contigs. If one contig is fully covered 
(95% threshold) by other contigs of the draft assembly, the contig is removed. 
blast+ is used in this script. If blast+ is not in the path, the path can
be added to the variable blast. 
'''
)
parser.add_argument("-d", "--draft", help="draft genome assembly (multi fasta)", required=True)
parser.add_argument("-o", "--output", help="Prefix of the output file", required=True)
parser.add_argument("-b", "--blastPath", help="Path to blast+ function, if not in path", type=str, default="")

# Read arguments from the command line
args = parser.parse_args()

draftPath=args.draft
outputPath=args.output
blast=args.blastPath
if blast != "" and not blast.endswith("/") :
	blast += "/"

print("\n\t--- REMOVING REDUNDANT CONTIGS ---\n")
print("Arguments detected:")
print(f"\t--draft:\t{draftPath}")
print(f"\t--output:\t{outputPath}\n")

# ===============
# Get Input files
# ===============
# Get draft contigs names and sequences
# Read draft Fasta
draftFasta = Fasta(draftPath)

blastResultsPath = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000)) + ".blastn"


# Run Blastn of draft against itself
start = time.time()
blastn(draftPath, draftPath, blast, blastResultsPath)
end = time.time()
print("Alignment done: ran in "+str(round(end-start))+"s")

# Get BED of the draft assemblies (BEDcoordinates without N nucleotides)
start = time.time()
draftBED = getBED(draftPath)
end = time.time()
print(f"Draft BED obtained: ran in {round(end-start)}s")

# Blast results analysis
# Blast results will be stored in the 2 level matrix analysisMatrix
# Level 1: Contig on which the contigs are aligned
# Level 2: BED of alignments covering contig1 by contig2
start = time.time()
x = len(draftFasta)
y = len(draftFasta)
analysisMatrix = []
for i in range(x):
	analysisMatrix.append([])
	for j in range(y):
		analysisMatrix[i].append([])

# Read blast out file
with open(blastResultsPath) as blastResultsFile:
	blastResults = csv.reader(blastResultsFile, delimiter="\t")
	for row in blastResults:
		contig1 = row[0]
		contig2 = row[1]
		if contig1 != contig2:
			index1 = draftFasta.getIndexFromID(contig1)
			index2 = draftFasta.getIndexFromID(contig2)
			startPos = int(row[6])
			endPos = int(row[7])+1
			coord2add = BEDcoordinates(id = contig1, start = startPos, end = endPos)
			analysisMatrix[index1][index2] += [coord2add]

os.remove(blastResultsPath)

# Convert every list of BEDcoordinates to BED objects
for i in range(len(draftFasta)):
	for j in range(len(draftFasta)):
		analysisMatrix[i][j] = BED(analysisMatrix[i][j])
end = time.time()
print("Blast read: ran in "+str(round(end-start))+"s")

# Round 1: Compute coverage of each contig by all other contigs
start = time.time()
alignments = []
overlap = []
for i in range(len(draftFasta)):
	# For each contig, the alignment list will contain the alignment BED of all other contigs on the contig
	contigs2align = list(range(len(draftFasta)))
	contigs2align.remove(i)
	BED2sum = [analysisMatrix[i][x] for x in contigs2align]
	alignments += [BED(BED2sum)]

	# Convert each BED to overlap
	draftContigBED = draftBED.getID(draftFasta.getID()[i]) # Get BED of contig
	overlap += [draftContigBED.overlapLen(alignments[i], percent = True)] # Get overlap percentage of alignment on contig

# While there is contigs covered on more than 95%, remove contig and recompute the whole process
contigsRemoved = []
contigsRemovedCoverage = []
while max(overlap) >= 95:
	contigsRemoved += [overlap.index(max(overlap))]
	contigsRemovedCoverage += [max(overlap)]

	alignments = [[] for i in range(len(draftFasta))]
	overlap = [0 for i in range(len(draftFasta))]
	# For each contig, the alignment list will contain the alignment BED of all other contigs on the contig
	for i in range(len(draftFasta)):
		if i not in contigsRemoved:
			# For each contig, the alignment list will contain the alignment BED of all other contigs on the contig
			contigs2align = list(range(len(draftFasta)))
			contigs2align.remove(i)
			for j in contigsRemoved:
				contigs2align.remove(j)
			BED2sum = [analysisMatrix[i][x] for x in contigs2align]
			alignments[i] = BED(BED2sum)

			# Convert each BED to overlap
			draftContigBED = draftBED.getID(draftFasta.getID()[i])
			overlap[i] = draftContigBED.overlapLen(alignments[i], percent = True)

end = time.time()
print("Coverage analysis completed: ran in "+str(round(end-start))+"s")

# Create Fasta objects with redundant and non redundant contigs
nonRedundantContigs = Fasta()
redundantContigs = Fasta()
for i in range(len(draftFasta)):
	seq = draftFasta.sequences[i]
	if i in contigsRemoved:
		redundantContigs += Fasta([seq])
	else:
		nonRedundantContigs += Fasta([seq])

# Write output to 2 files: PREFIX.NoRedundantContigs.fasta and PREFIX.removedContigs.fasta
nonRedundantContigs.toFile(outputPath+".NoRedundantContigs.fasta")

if len(contigsRemoved) > 0 :
	redundantContigs.toFile(outputPath+".removedContigs.fasta")

if len(contigsRemoved) > 0 :
	print("\nContigs removed: ")
	print("\tContig\t% covered")
	for seq in redundantContigs:
		print(f"\t{seq.id}\t{contigsRemovedCoverage[redundantContigs.sequences.index(seq)]}")
else:
	print("\nNo contig removed")
