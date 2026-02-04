#!/bin/env/python3

from datetime import datetime as dt
import pandas as pd
import time
import argparse

## function to convert datetime object to decimal
def toYearFraction(date):
    def sinceEpoch(date): # returns seconds since epoch
        return time.mktime(date.timetuple())
    s = sinceEpoch

    year = date.year
    startOfThisYear = dt(year=year, month=1, day=1)
    startOfNextYear = dt(year=year+1, month=1, day=1)

    yearElapsed = s(date) - s(startOfThisYear)
    yearDuration = s(startOfNextYear) - s(startOfThisYear)
    fraction = yearElapsed/yearDuration

    return date.year + fraction

## parse command-line options
parser = argparse.ArgumentParser(description='Sort metadata by date and write to outfile')
parser.add_argument('-i', '--infile', metavar='INPUT FILE', action='store', type=str, required=True, help='metadata of genome samples in .tsv format (with date attribute for sampling date')
parser.add_argument('-o', '--outfile', metavar='OUTPUT FILE', action='store', type=str, required=True, help='output file to be written to')
parser.add_argument('-l', '--latest', action='store_true', default=False, help='output last sample date in decimal only')
args = parser.parse_args()

## 
meta = pd.read_csv(args.infile, sep='\t')
meta['date_dec'] = meta.date.apply(lambda x: toYearFraction(dt.strptime(x, '%Y-%m-%d')))
meta = meta.sort_values(by=['date_dec'], ascending=False)

## print last sample date in decimal if specified
if args.latest:
    print(meta.date_dec.max())
else: ## otherwise write full metadata to outfile
    meta.to_csv(args.outfile, sep='\t', index=False)