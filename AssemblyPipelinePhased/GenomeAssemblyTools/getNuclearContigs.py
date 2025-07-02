#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/01/04
# version ='2.0'
# ---------------------------------------------------------------------------
'''
This script retrieve nuclear contigs and discard small contigs. A nuclear 
contigs is defined as a contig that aligns on the reference nuclear genome
over more than X consecutive kb with blast. 

It takes as input :
	-d --draft: a draft genome assembly to reorder (multi fasta)
	-db --database: the blast database of nuclear chromosomes (done with makeblastdb)
	-r --reference: reference of the nuclear genome (multi fasta)
	-m --minLength: minimum length of a contig (kb)
	-a --minAlignedLength: minimum length of the alignment of the contig on the ref (kb)
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
def blastnDB(query, database, blastPath, out):
	blastCommand = blastPath + "blastn -query " + query + " -db " + database + " -dust no -soft_masking false -max_target_seqs 1 -max_hsps 1 -num_threads 4 -outfmt '6 qseqid qlen length' -out " + out
	print(blastCommand)
	os.system(blastCommand)
def blastnRef(query, reference, blastPath, out):
	blastCommand = blastPath + "blastn -query " + query + " -subject " + reference + " -dust no -soft_masking false -max_target_seqs 1 -max_hsps 1 -num_threads 4 -outfmt '6 qseqid qlen length' -out " + out
	print(blastCommand)
	os.system(blastCommand)
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(
	description = '''
This script retrieve nuclear contigs and discard small contigs. A nuclear 
contigs is defined as a contig that aligns on the reference nuclear genome
over more than X consecutive kb with blast. 
blast+ is used in this script. If blast+ is not in the path, the path can
be added to the variable blast. 
'''
)
parser.add_argument("-d", "--draft", help="draft genome assembly (multi fasta)", required=True)
parser.add_argument("-db", "--database", help="the blast database of nuclear chromosomes (done with makeblastdb). Either reference or Database must be given", default = "")
parser.add_argument("-r", "--reference", help="referance of the nuclear genome (multifasta). Either reference or Database must be given", default = "")
parser.add_argument("-m", "--minLength", help="minimum length of a contig (kb)", type=int, default=1)
parser.add_argument("-a", "--minAlignedLength", help="minimum alignment length on the nuclear genome (kb)", type=float, default=2)
parser.add_argument("-o", "--output", help="Prefix of the output file", required=True)
parser.add_argument("-b", "--blastPath", help="Path to blast+ function, if not in path", type=str, default="")

# Read arguments from the command line
args = parser.parse_args()

draftPath=args.draft
dbPath=args.database
refPath=args.reference
minLength=args.minLength
minAlignedLength=args.minAlignedLength
outputPath=args.output
blast=args.blastPath
if blast != "" and not blast.endswith("/") :
	blast += "/"
if refPath == "" and dbPath == "":
	raise ValueError("You must either specify the path to a blastn database or to a reference in fasta format")
elif refPath != "" and dbPath != "":
	raise ValueError("You cannot specify both the path to a blastn database and to a reference in fasta format")



print("\n\t--- KEEPING NUCLEAR CONTIGS ---\n")
print("Arguments detected:")
print(f"\t--draft:\t{draftPath}")
if dbPath != "":
	print(f"\t--database:\t{dbPath}")
else:
	print(f"\t--reference:\t{refPath}")	
print(f"\t--minLength:\t{minLength}kb")
print(f"\t--minAlignedLength:\t{minAlignedLength}kb")
print(f"\t--output:\t{outputPath}\n")

# ===============
# Get Input files
# ===============
# Read draft fasta
draftFasta = Fasta(draftPath)

blastResultsPath = datetime.now().strftime("%Y_%m_%d_%H_%M_%S_%m_%f") + "_" + str(randint(0, 10000)) + ".blastn"

# Run Blastn of draft against ref
start = time.time()
if dbPath != "":
	blastnDB(draftPath, dbPath, blast, blastResultsPath)
else:
	blastnRef(draftPath, refPath, blast, blastResultsPath)	
end = time.time()
print("Alignment done: ran in "+str(round(end-start))+"s")

# Blast result analysis
# Keep contigs with length >= minLength and that aligned on nuclear ref genome
# for >= minAlignedLength consecutive kb
contigs2keep = []
blastResultsFile = open(blastResultsPath)
blastResults = csv.reader(blastResultsFile, delimiter="\t")
for row in blastResults:
	contig = row[0]
	contigLength = int(row[1])
	alignmentLength = int(row[2])
	if contigLength >= minLength * 1000 and alignmentLength >= minAlignedLength * 1000:
		contigs2keep += [contig]
blastResultsFile.close()
os.remove(blastResultsPath)

# Split contigs between nuclear and non nuclear
nuclearFasta = Fasta()
nonNuclearFasta = Fasta()
for seq in draftFasta:
	if seq.id in contigs2keep:
		nuclearFasta += Fasta([seq])
	else:
		nonNuclearFasta += Fasta([seq])

# write output file
nuclearFasta.toFile(outputPath+".Nuclear.fasta")

if len(nonNuclearFasta) > 0 :
	nonNuclearFasta.toFile(outputPath+".NonNuclear.fasta")

print(f'\n{len(nuclearFasta)}/{len(draftFasta)} contigs kept as NUCLEAR CONTIGS\n')
