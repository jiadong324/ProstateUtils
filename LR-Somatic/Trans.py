#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/5/17

'''


import math

import matplotlib.pylab as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import intervaltree
import seaborn as sns
from matplotlib.patches import Patch
import upsetplot

from Helpers.Constants import *
from Helpers.Functions import *

sns.set_theme(style="ticks", font="Arial", font_scale=1.0)

plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2


def filter_trans_bp(workdir):
    sds = pysam.Tabixfile(REMOTESD, 'r')
    ref_file = pysam.FastaFile(REMOTEREF)

    for sample in ['E']:
        normal_bps = {}
        for line in open(f'{workdir}/{sample}-B.s3.old_format.bed', 'r'):
            entries = line.strip().split('\t')
            trans_id, chrom, start, end = entries[0], entries[3], int(entries[4]) - 100, int(entries[4]) + 100
            if chrom in normal_bps:
                normal_bps[chrom][start:end] = (start, end)
            else:
                normal_bps[chrom] = IntervalTree()
                normal_bps[chrom][start:end] = (start, end)

        somatic_trans = open(f'{workdir}/{sample}.somatic.s3.old_format.bed', 'w')
        normal_id = []
        somatic_bps = []
        counter_id = set()

        for line in open(f'{workdir}/{sample}-T.s3.old_format.bed', 'r'):
            entries = line.strip().split('\t')
            trans_id, chrom, start, end = entries[0], entries[3], int(entries[4]) - 100, int(entries[4]) + 100

            if chrom not in VALID_CHROMS:
                continue

            if normal_bps[chrom].overlap(start, end):
                normal_id.append(trans_id)
                continue

            if contains_gaps(chrom, start, end, ref_file):
                normal_id.append(trans_id)
                continue

            for sd in sds.fetch(chrom, start, end):
                entries = sd.strip().split('\t')
                sd_start, sd_end, sd_mate_coord = int(entries[1]), int(entries[2]), entries[3]
                overlap_size = get_overlaps(start, end, sd_start, sd_end)
                if overlap_size > 0:
                    normal_id.append(trans_id)
                    continue

            somatic_bps.append([trans_id, line.strip()])

        for trans_info in somatic_bps:
            if trans_info[0] in normal_id:
                continue
            counter_id.add(trans_info[0])
            print(trans_info[1], file=somatic_trans)
        print(f'{sample}: {len(counter_id)}')

def contains_gaps(chrom, start, end, ref):
    seq = ref.fetch(chrom, start - 2000, end + 2000)
    in_gap = False
    for base in seq:
        if base == 'N':
            in_gap = True
            break
    return in_gap

def pieplot_trans_bp(workdir):

    bps_count = [0, 0, 0, 0]

    for i, sample in enumerate(SAMPLES):
        filtered_tras = pd.read_csv(f'{workdir}/somatic_trans/{sample}-T.trans.s2.ex_B.ex_in2.ex_Re.bed', sep='\t', header=[0])
        trans_dict = {}
        for idx, row in filtered_tras.iterrows():
            id, chrom, start, end, tratype = row['TransId'], row['Chr'], row['Start'], row['End'], row['TransType']
            if id in trans_dict:
                trans_dict[id].append((chrom, start, end))
            else:
                trans_dict[id] = [(chrom, start, end)]

        # this_sample_count = [0, 0, 0]

        for id, links in trans_dict.items():
            bps = len(links)
            if bps == 2:
                bps_count[0] += 1
            elif bps == 3:
                bps_count[1] += 1
            elif bps == 4:
                bps_count[2] += 1
            elif bps > 4:
                bps_count[3] += 1

        # bps_count.append(this_sample_count)
    explode = [0, 0, 0, 0.2]

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.pie(bps_count, labels=['2', '3', '4', '>4'], labeldistance=1.15, autopct=make_autopct(bps_count),
           wedgeprops={'linewidth': 1, 'edgecolor': 'white', 'width':0.5}, explode=explode)

    # bottoms = []
    # for j, sample in enumerate(SAMPLES):
    #     this_sample_data = bps_count[j]
    #
    #     bottoms.append(this_sample_data)
    #
    #     if j == 0:
    #         ax.bar([0, 1, 2], this_sample_data, label=sample, width=0.6, edgecolor='w', color=SAMPLECOLOR[sample])
    #     else:
    #         bottom_sum = [sum(x) for x in zip(*bottoms[0: j])]
    #         ax.bar([0, 1, 2], this_sample_data, bottom=bottom_sum, label=sample, width=0.6, color=SAMPLECOLOR[sample])
    #
    # ax.set_xticks([0, 1, 2])
    # ax.set_xticklabels(['3', '4', '>4'], fontsize=12)
    # ax.legend()
    # ax.set_ylabel('Number of links', fontsize=13)
    #
    # ax.spines['top'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    #
    # ax.set_yticks(np.linspace(0, 120, 5))
    # ax.set_yticklabels([int(val) for val in np.linspace(0, 120, 5)], fontsize=12)
    #
    # print(bps_count)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/somatic_trans_bp_count.pdf')
    plt.show()

def main():

    filter_trans_bp('/data/home/jdlin/Prostate/somatic_trans/ngmlr_s3_v0524')
    # pieplot_trans_bp('/data/home/jdlin/Prostate')

if __name__ == '__main__':
    main()

