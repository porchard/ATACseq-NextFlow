#!/usr/bin/env python
# coding: utf-8

import sys
import pybedtools as bt
import argparse

parser = argparse.ArgumentParser(description='Convert a set of summits to a set of peaks as described in methods of paper "Epigenomic State Transitions Characterize Tumor Progression in Mouse Lung Adenocarcinoma"', add_help=True)
parser.add_argument('--extension', type=int, default=150, help='Extend summits by this amount on either side (default: 150)')
parser.add_argument('summits_bed', type=str, help='Bed file of summits. Chrom, start, end, name (must be unique), -log10(q). As output from MACS2 callpeak --call-summits')
parser.add_argument('chrom_sizes', type=str, help='Chrom size file')
args = parser.parse_args()

SUMMITS = args.summits_bed
EXTENSION = args.extension
CHROM_SIZES = args.chrom_sizes


if not EXTENSION >= 0:
    raise ValueError('--extension must be a value greater than 0')


# extend summits --> peaks
peaks = bt.BedTool(SUMMITS).sort().slop(b=EXTENSION, g=CHROM_SIZES).saveas()

# find overlapping peaks
overlaps = peaks.intersect(peaks, wa=True, wb=True).to_dataframe()
overlaps.columns = ['chrom_1', 'start_1', 'end_1', 'name_1', 'score_1', 'chrom_2', 'start_2', 'end_2', 'name_2', 'score_2']
overlaps = overlaps[overlaps.name_1 != overlaps.name_2]

if len(overlaps) > 0:
    overlapping = {i: df.name_2.to_list() for i, df in overlaps.groupby('name_1')} # peak --> [overlapping_peak_1, overlapping_peak_2, ...]
else:
    overlapping = {}


# drop less-significant overlapping peaks
peaks = peaks.to_dataframe().sort_values('score', ascending=False)
drop_peaks = set() # drop these peaks in the end, as they overlap with peaks from more significant summits

for peak in peaks.name:
    if peak in drop_peaks: # this peak overlaps a more significant peak -- so ignore it
        continue
    if peak not in overlapping: # there are no peaks overlapping this one, so keep it
        continue
    else:
        # each of the overlapping peaks need to be marked as DROP
        for overlapping_peak in overlapping[peak]:
            drop_peaks.add(overlapping_peak)

sys.stderr.write('Dropping {:,} of {:,} peaks\n'.format(len(drop_peaks), len(peaks)))
peaks[~peaks.name.isin(drop_peaks)].to_csv(sys.stdout, sep='\t', index=False, header=None)