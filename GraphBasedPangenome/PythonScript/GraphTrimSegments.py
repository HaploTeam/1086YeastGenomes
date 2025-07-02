#!/usr/bin/env python
# -*- coding: utf8 -*-
#----------------------------------------------------------------------------
# Created By  : vloegler
# Created Date: 2024/02/28
# version ='1.0'
# ---------------------------------------------------------------------------
'''
This script remove short segments from a graph at the gfa format.
'''
# ---------------------------------------------------------------------------
import pandas as pd
import argparse
import os
import re
import itertools
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
parser.add_argument("-g", "--gfa", help="Path of the graph (gfa format)", required=True)
parser.add_argument("-s", "--min_size", help="Minimum segment size (bp)", required=True, type = int)

# Read arguments from the command line
args = parser.parse_args()

graphPath = os.path.abspath(args.gfa)
min_size = args.min_size

#graphPath = '/ccc/scratch/cont007/fg0006/loeglerv/GraphPangenome/04-analysis/12_GraphUniqueSequences/SubGraph.s64709.c50.gfa'

# Read the graph
segments = []
links = []
with open(graphPath, 'r') as graph:
    for line in graph:
        line = line.strip()
        if line.startswith('S'):
            segments += [pd.DataFrame(data = {'ID': [line.split('\t')[1]], 'Seq': [line.split('\t')[2]]})]
        elif line.startswith('L'):
            links += [pd.DataFrame(data = {'ID1': [line.split('\t')[1]], 'Direction1': [line.split('\t')[2]], 'ID2': [line.split('\t')[3]], 'Direction2': [line.split('\t')[4]]})]

segments = pd.concat(segments).reset_index(drop=True)
segments['Len'] = segments['Seq'].apply(len)
links = pd.concat(links).reset_index(drop=True)

# Iterate on segments shorter than minimum size
for i in segments.loc[segments['Len'] < min_size,:].index:
    print(i)
    # Get segment
    seg = segments.loc[i, 'ID']
    # Get all segments linked
    upstream_linked_seg = [(links.loc[i,'ID1'], links.loc[i,'Direction1']) for i in links.loc[(links['ID2'] == seg), :].index]
    downstream_linked_seg = [(links.loc[i,'ID2'], links.loc[i,'Direction2']) for i in links.loc[(links['ID1'] == seg), :].index]
    # Remove segment
    segments = segments.loc[segments['ID'] != seg,:]
    # Remove links including segment
    links = links.loc[(links['ID1'] != seg) & (links['ID2'] != seg),:]
    print(f'{len(upstream_linked_seg)} upstreamSeg and {len(downstream_linked_seg)} down')
    if len(upstream_linked_seg) > 0 and len(downstream_linked_seg) > 0:
        # Add new links
        up = list(itertools.chain.from_iterable([[x] * len(downstream_linked_seg) for x in upstream_linked_seg]))
        down = downstream_linked_seg * len(upstream_linked_seg)
        links = pd.concat([links, pd.DataFrame(data = {'ID1': [x[0] for x in up], 'Direction1': [x[1] for x in up], 'ID2': [x[0] for x in down], 'Direction2': [x[1] for x in down]})]).reset_index(drop = True).drop_duplicates()

# Write the new gfa file
with open(f'{re.sub(".gfa", "", graphPath)}.Trimmed{min_size}bp.gfa', 'w') as out:
    # Write header
    out.write('H\tVN:Z:1.1\n')

# Write segments
segments['name'] = 'S'
segments.loc[:,['name', 'ID', 'Seq']].to_csv(f'{re.sub(".gfa", "", graphPath)}.Trimmed{min_size}bp.gfa', index = False, header = False, mode = 'a', sep = '\t')

# Write links
links['name'] = 'L'
links['end'] = '*'
links.loc[:,['name', 'ID1', 'Direction1', 'ID2', 'Direction2', 'end']].to_csv(f'{re.sub(".gfa", "", graphPath)}.Trimmed{min_size}bp.gfa', index = False, header = False, mode='a', sep = '\t')






