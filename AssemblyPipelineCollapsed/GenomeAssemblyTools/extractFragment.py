#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2022/09/15
# version ='2.1'
# ---------------------------------------------------------------------------
'''
This script extract a sequence from a fasta file

It takes as input :
	-f --fasta: fasta file
	-F --file: Text file containing IDs to extract, one ID per line
	-id --seqid: ID of the sequence, or multiple IDs separated by space
	-p --prefix: Prefix that the sequence ID have to match to be retained (several allowed seperated by space)
	-s --start: Start position (inclusive, 1-based)
	-e --end: End position (inclusive, 1-based)
	-o --output: output file
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
parser = argparse.ArgumentParser(description = "This script extracts a sequence from a fasta file. It can also extract a fragment of a sequence, retrieve all sequences whose ID starts with a specified prefix, or retrieve the opposite of a query. ")
parser.add_argument("-f", "--fasta", help="fasta file", required=True)
parser.add_argument("-F", "--file", help="Text file containing IDs to extract, one ID per line", default = "")
parser.add_argument("-o", "--output", help="output file", default = "")
parser.add_argument("-id", "--seqid", help="ID of the sequence, or multiple IDs separated by space", nargs='+', default = [])
parser.add_argument("-p", "--prefix", help="Prefix that the sequence ID have to match to be retained (several allowed seperated by space)", nargs='+', default = [])
parser.add_argument("-s", "--start", help="Start position (inclusive, 1-based)", type = int, default = 0)
parser.add_argument("-e", "--end", help="End position (inclusive, 1-based)", type = int, default = 0)
parser.add_argument("-v", "--invert", help="Output the inverse of the query", action='store_true')

# Read arguments from the command line
args = parser.parse_args()

fastaPath = args.fasta
IDpath = args.file
outputPath=args.output
seqidList=args.seqid
if IDpath != "":
	# If a file of IDs is given, add IDs to seqidlist
	with open(IDpath, "r") as IDs:
		seqidList += [line.strip() for line in IDs]
prefixList=args.prefix
start=int(args.start)
end=int(args.end)

# ===============
# Get Input files
# ===============

# Read fasta file
fasta = Fasta(fastaPath)

# Open output file if necessary
if outputPath != "":
	out = open(outputPath, "w")

# If coordinates are present, check proper coordinates 
if start != 0 or end != 0:
	if start < 1 or end < 1 or end < start:
		raise ValueError("Bad coordinates. Start and End must be 1-based inclusive and start <= end.")


# Classic mode ===============================================================

if not args.invert:
	# Check if all sequence ID or prefix are present in the fasta to extract
	for seqid in seqidList:
		if seqid not in fasta.getID():
			raise ValueError("Sequence ID " + seqid + " is not in fasta file.")
	for prefix in prefixList:
		if not any(item.startswith(prefix) for item in fasta.getID()):
			raise ValueError("Prefix " + prefix + " is not in fasta file.")

	for ID in fasta.getID():
		# Check if ID is in the seqid list 
		# 	or starts with any prefix in prefixList
		# 	or no ID nor prefix was given
		# 	and start of the sequence is smaller than sequence length
		if ((ID in seqidList) or any(ID.startswith(item) for item in prefixList) or len(seqidList)+len(prefixList) == 0) and start <= len(fasta.getSeqFromID(ID)):

			# Retrieve sequence corresponding to the ID
			targetSeq = fasta.getSeqFromID(ID)

			# If coordinates are present, extract the fragment
			if start != 0 or end != 0:
				# if end is higher than fragment length, scale down end
				if end > len(targetSeq):
					endPos = len(targetSeq)
				else:
					endPos = end
				# Extract fragment
				extractedSeq = Sequence(f'{ID}_{start}..{endPos}', targetSeq[(start-1):endPos])
			else: # If no coordinate given
				extractedSeq = Sequence(ID, targetSeq)
			

			# Write the sequence to the ouptupt (stdout or file)
			if outputPath != "":
				out.write(extractedSeq.__str__()+"\n")
			else:
				sys.stdout.write(extractedSeq.__str__()+"\n")


# Invert mode ===============================================================

if args.invert:

	for ID in fasta.getID():
		# Check if ID is not in the seqid list 
		# 	and doesn't start with any prefix in prefixList
		# 	and start-end does not contain the whole sequence
		if (ID not in seqidList) and not any(ID.startswith(item) for item in prefixList) and not (start == 1 and end >= len(fasta.getSeqFromID(ID))):

			# Retrieve sequence corresponding to the ID
			targetSeq = fasta.getSeqFromID(ID)
			
			# If coordinates are present, take the inverse and extract the fragment
			if start != 0 or end != 0:
				# if start begin at 1, extract seq from end to end of chr
				if start == 1:
					nFragment = 1 # Nb of remaining fragments
					startPos = [end+1]
					endPos = [len(targetSeq)]
				# if end finish at the end of chr, extract seq from the begining to start
				elif end >= len(targetSeq):
					nFragment = 1 # Nb of remaining fragments
					startPos = [1]
					endPos = [start-1]
				# if start and end are comprised in the sequence, split the sequence in 2
				else:
					nFragment = 2 # Nb of remaining fragments
					startPos = [1, end+1]
					endPos = [start-1, len(targetSeq)]

				# Extract fragment
				extractedSeq = []
				for i in range(nFragment):
					extractedSeq += [Sequence(f'{ID}_{startPos[i]}..{endPos[i]}', targetSeq[(startPos[i]-1):endPos[i]])]
				extractedSeq = Fasta(extractedSeq) # Convert to Fasta object
			else: # If no coordinate given
				extractedSeq = Fasta([Sequence(ID, targetSeq)])

			# Write the sequence to the ouptupt (stdout or file)
			if outputPath != "":
				out.write(extractedSeq.__str__()+"\n")
			else:
				sys.stdout.write(extractedSeq.__str__()+"\n")

if outputPath != "":
	out.close()

