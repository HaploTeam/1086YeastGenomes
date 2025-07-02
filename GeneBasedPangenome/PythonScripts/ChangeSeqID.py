#!/usr/bin/env python
# -*- coding: utf8 -*-
# This script change IDs of the sequences of a fasta file. 
# It takes a prefix and each sequence is named Prefix_$i, 
# with i increasing from 1.

import argparse
import re
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser(description = """
This script change IDs of the sequences of a fasta file. 
It takes a prefix and each sequence is named Prefix_$i, 
with i increasing from 1.
""")
parser.add_argument("-f", "--fasta", help="fasta file", required=True)
parser.add_argument("-p", "--prefix", help="prefix", required=True)
parser.add_argument("-o", "--output", help="file output", required=True)

# Read arguments from the command line
args = parser.parse_args()

outputPath = args.output
prefix = args.prefix
fastaPath = args.fasta

c = 1

out = open(outputPath, "w")
with open(fastaPath, "r") as fasta:
    for line in fasta:
        if line.startswith('>'):
            newline = re.sub('>', f'>{prefix}_{c} ', line)
            c += 1
            out.write(newline)
            #print(newline)
        else:
            out.write(line)
            #print(line)
out.close()
