#!/bin/env/python3

import baltic as bt
import argparse

## parse command-line options
parser = argparse.ArgumentParser(description='Annotate and count number of state changes from results of TreeTime mugration')
parser.add_argument('-i', '--infile', metavar='INPUT FILE', action='store', type=str, required=True, help='unprocessed output file TreeTime mugration (with country as attribute')
parser.add_argument('-o', '--outfile', metavar='OUTPUT FILE', action='store', type=str, required=True, help='output file to be written to')
parser.add_argument('-t', '--time', metavar='LATEST DATE', action='store', type=float, required=True, help='time of last sampled tip')

args = parser.parse_args()

tree = bt.loadNewick(args.infile, absoluteTime=False)
tree.setAbsoluteTime(args.time)

tree.traverse_tree() ## required to set heights
tree.treeStats() ## report stats about tree

changes = 0
times = []
origins = []
destinations = []
for k in tree.Objects: ## iterate over a flat list of branches
    
    # Assign node UNKNOWN country if not give
    if 'country' in k.traits:
        country = k.traits['country']
    else:
        country = 'UNKNOWN'
        k.traits['country'] = country 
    
    # Find parent country if given
    if k.parent.traits:
        parent_country = k.parent.traits['country']
    else:
        parent_country = 'UNKNOWN'
        
    if country != parent_country:
        changes += 1
        times = times + [k.absoluteTime]
        origins = origins + [parent_country]
        destinations = destinations + [country]

with open(args.outfile, 'w+') as outfile:
    outfile.write('EventTime\tOrigin\tDestination\n')
    outfile.write('\n'.join(['%f\t%s\t%s' % (time, origins[i], destinations[i]) for i, time in enumerate(times)]))