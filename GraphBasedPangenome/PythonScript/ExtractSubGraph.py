#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2024/02/28
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script extract a snarl from a pangenome graph (gbz) using the variant ID
found in the graph VCF. 
It uses bcftools and vg. 
'''
# ---------------------------------------------------------------------------
import pandas as pd
import argparse
import re
import os
import subprocess
import logging
logging.basicConfig(level=logging.NOTSET)
# ---------------------------------------------------------------------------
# Definitions
# ---------------------------------------------------------------------------

# =============
# Get arguments
# =============

# Initiate the parser
parser = argparse.ArgumentParser()
parser.add_argument("-v", "--vcf", help="Path of the graph VCF", required=True)
parser.add_argument("-i", "--id", help="ID of the variant to extract", required=True)

parser.add_argument("-g", "--graph", help="Path of the graph gbz file", required=True)
parser.add_argument("-o", "--output", help="Output file (gfa format)", required=True)

# Read arguments from the command line
args = parser.parse_args()

vcfPath=os.path.abspath(args.vcf)
id=args.id
graphPath=os.path.abspath(args.graph)
outputPath=os.path.abspath(args.output)

#graphPath = '/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis/06_ReferenceGraph/SacePangenomeGraph.500Haplotypes.giraffe.gbz'
#vcfPath = '/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis/10_PopulationVCF/GraphGenotyping.2874Samples.chromosome7.vcf.gz'
#id = '>11031218>11031250'

# Extract the INFO/AT tag of the variant
logging.info('Extracting variants from the VCF file')
cmd1 = ['bcftools', 'view', 
        vcfPath, 
        '-i', f'ID==\'{id}\'']
cmd2 = ['bcftools', 'query', 
        '-f', '"%INFO/AT\n"']
#ps = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
#paths = subprocess.check_output(cmd2, stdin=ps.stdout)
#paths = re.sub('"', '', paths.decode("utf-8")).strip()

cmd1 = ['zcat', vcfPath]
cmd2 = ['sed', f'/\t{id}\t/q']
cmd3 = ['tail', '-1']
cmd4 = ['cut', '-f', '8']
c1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
c2 = subprocess.Popen(cmd2, stdin=c1.stdout, stdout=subprocess.PIPE)
c3 = subprocess.Popen(cmd3, stdin=c2.stdout, stdout=subprocess.PIPE)
paths = subprocess.check_output(cmd4, stdin=c3.stdout)
paths = re.sub('"', '', paths.decode("utf-8")).strip()
paths = [x for x in paths.split(';') if x.startswith('AT=')][0]
paths = paths.split('=')[1]

#zcat 10_PopulationVCF/GraphGenotyping.2874Samples.chromosome2.vcf.gz | sed '/\t>6706081>6706092\t/q' | tail -1 | cut -f 8

#paths='>6706081>6706082>6706084>6706085>6706088>6706089>6706092,>6706081>6706082>6706084>6706085>6706088>6706091>6706092,>6706081>6706082>6706084>6706087>6706088>6706089>6706092,>6706081>6706082>6706084>6706087>6706088>6706091>6706092,>6706081>6706082>6706084>6706087>6706090>6706091>6706092,>6706081>6706083>6706084>6706087>6706088>6706089>6706092,>6706081>6706083>6706084>6706087>6706088>6706091>6706092,>6706081>6706083>6706084>6706085>6706088>6706089>6706092,>6706081>6706083>6706084>6706087>6706090>6706091>6706092,>6706081>6706083>6706084>6706085>6706088>6706091>6706092,>6706081>6706083>6706086>6706087>6706088>6706089>6706092,>6706081>6706083>6706086>6706087>6706090>6706091>6706092'

# Write header ================================================
with open(outputPath, 'w') as out:
        out.write(f'H\tVN:Z:1.1\tRS:Z:S288c\tVARIANT:{id}\n')

# Get segments ================================================
logging.info('Extracting segments from graph')
# Get the unique list of segments
segments = re.sub('>', ',', re.sub('<', ',', paths)).split(',')
segments = [x for x in segments if x != '']
segments = sorted(set(segments))
# Retrieve segments from graph
with open('segments.tmp', 'w') as tmp:
    for s in segments:
        tmp.write(f'[^0-9]{s}[^0-9]\n')
cmd1 = ['vg', 'view', graphPath]
cmd2 = ['sed', "/^W/q"]
cmd3 = ['grep', '^S']
cmd4 = ['grep', '-f', 'segments.tmp']
c1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
c2 = subprocess.Popen(cmd2, stdin=c1.stdout, stdout=subprocess.PIPE)
c3 = subprocess.Popen(cmd3, stdin=c2.stdout, stdout=subprocess.PIPE)
segmentLines = subprocess.check_output(cmd4, stdin=c3.stdout)
segmentLines = re.sub('"', '', segmentLines.decode("utf-8"))
# Write to file
with open(outputPath, 'a') as out:
        out.write(segmentLines)

# Get links ========================================================
logging.info('Extracting links from variants')
links = []
for variant in paths.split(','):
    variant = re.sub('>', ',>', re.sub('<', ',<', variant))
    variant = [(x[1:], '+') if x[0] == '>' else (x[1:], '-') for x in variant.split(',') if x != '']
    for i in range(len(variant) - 1):
        links += [pd.DataFrame(data = {'RecordType': ['L'], 'From': [variant[i][0]], 'FromOrient': [variant[i][1]], 'To': [variant[i+1][0]], 'ToOrient': [variant[i+1][1]], 'Overlap': ['*'], })]

links = pd.concat(links).drop_duplicates()

# Write to file
links.to_csv(outputPath, sep = '\t', index = False, header = False, mode = 'a')
logging.info(f'Graph successfully extracted to {outputPath}')

