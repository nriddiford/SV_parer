#!/usr/bin/env python

import sys
import fnmatch
import os
from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f",
                  "--SV_calls",
                  dest="SV_calls",
                  help="SV calls file"
                  )

options, args = parser.parse_args()


if not options.SV_calls:
    parser.error('No input file provided')

variants_file = options.SV_calls

base_name = (os.path.splitext(variants_file)[0])
out_base = base_name.split('_')[0]

pattern = out_base + '*.cnv'

big_window='/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data'
samll_window='/Users/Nick_curie/Desktop/script_test/cnvPlotteR/data/w500'

def select_cnvfile(sv_size, bp1, bp2):
    if float(sv_size) <= 5000:
        files = os.listdir(samll_window)
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = "data/w500/" + name
                return[name, start, stop, tick]
    else:
        files = os.listdir(big_window)
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                name = "data/" + name
                return[name, start, stop, tick]

def print_del_depth(SV_calls_file):
    i=1
    with open(SV_calls_file, 'U') as calls:
        for l in calls:
            parts = l.rstrip().split('\t')

            if i == 1 and parts[0] == 'event':
                continue
            elif parts[9] != 'somatic_tumour':
                continue

            event, source, type, chrom1, bp1, chrom2, bp2 = parts[0:7]
            bp1 = int(bp1)
            bp2 = int(bp2)
            size = parts[11]

            if chrom1 == chrom2:
                fields = select_cnvfile(size, bp1, bp2)
                cnv_file = fields[0]
                start = fields[1]
                end = fields[2]
                tick = fields[3]
                print("SV event: %s, type: %s, size: %s") % (event, type, size)
                # print("regionPlot(cnv_file=\"%s\", from=%s, to=%s, bp1=%s,bp2=%s,chrom=\"%s\", tick=%s, title=\"%sKb %s on %s\")") % (cnv_file, start, end, bp1, bp2, chrom1, tick, size, type, chrom1)

print_del_depth(variants_file)
