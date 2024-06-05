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
import pysam
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

def plot_ins_info(workdir):

    ins_size = []

    for sample in SAMPLES:
        vcf_file = f'{workdir}/{sample}/sniffles/GRCh38/{sample}.somatic.aligners.jasmine.vcf'
        for line in open(vcf_file, 'r'):
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_dict = {}
            for token in entries[7].split(";"):
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            svlen = abs(int(info_dict['SVLEN']))
            if svlen > 100000:
                continue

            if 'SVTYPE' in info_dict:
                svtype = info_dict['SVTYPE']
            else:
                svtype = entries[2].split('-')[-2]

            if svtype == 'INS':
                ins_size.append([svlen, sample])

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    df_ins_size = pd.DataFrame(ins_size, columns=['svlen', 'sample'])
    sns.histplot(data=df_ins_size, x='svlen', log_scale=True, hue='sample', alpha=1, bins=50,
                 hue_order=SAMPLES, palette=[SAMPLECOLOR[ele] for ele in SAMPLES], ax=ax)

    ax.set_yticks(np.linspace(0, 60, 4))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 60, 4)], fontsize=12)
    ax.set_ylabel('# of INS', fontsize=13)

    ax.set_xlabel('SV length (bp)', fontsize=13)
    ax.set_xticks([50, 100, 300, 1000, 6000, 10000])
    ax.set_xticklabels(['50', '100', '300', '1,000', '6,000', '10,000'], fontsize=12, rotation=60)

    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=0)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/samples_somatic_ins_length.pdf')

    rm_dict = {}
    for line in open(f'{workdir}/sniffles_somatic/RM_INS_results/samples.somatic.ins.rm.txt', 'r'):
        entries = line.strip().split('\t')
        rm_type = entries[6]

        if rm_type in rm_dict:
            rm_dict[rm_type].append([entries[0], entries[7]])
        else:
            rm_dict[rm_type] = [[entries[0], entries[7]]]

    sorted_rm_counts = sorted(rm_dict.items(), key=lambda x:len(x[1]), reverse=True)
    labels = []
    counts = []
    for ele in sorted_rm_counts:
        if len(ele[1]) > 20:
            labels.append(ele[0])
            counts.append(len(ele[1]))

    explode = [0, 0, 0, 0, 0.2, 0, 0.2]
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    ax.pie(counts, labels=labels, labeldistance=1.15, autopct=make_autopct(counts), explode=explode, wedgeprops={'linewidth': 1, 'edgecolor': 'white'})
    plt.tight_layout()
    fig.savefig(f'{REMOTEFIG}/samples_somatic_ins_rmtype.pdf')

    subtype_dict = {}
    for key in ['LINE/L1', 'SINE/Alu']:
        subtype_dict[key] = {}
        for ele in rm_dict[key]:
            sample_idx = SAMPLES.index(ele[0])
            if ele[1] in subtype_dict[key]:
                subtype_dict[key][ele[1]][sample_idx] += 1
            else:
                subtype_dict[key][ele[1]] = [0, 0, 0, 0, 0]
                subtype_dict[key][ele[1]][sample_idx] += 1

    fig, axes = plt.subplots(1, 2, sharex='col', sharey='row', figsize=(9, 4))
    fig_idx = 0

    for key, subtypes in subtype_dict.items():
        ax = axes[fig_idx]
        bottoms = []
        sorted_subtypes = sorted(subtypes.items(), key=lambda x: sum(x[1]), reverse=True)
        xlabels = [ele[0] for ele in sorted_subtypes]
        xticks = np.arange(len(xlabels))
        for j, sample in enumerate(SAMPLES):
            this_sample_num = []
            for info in sorted_subtypes:
                this_sample_num.append(info[1][j])

            bottoms.append(this_sample_num)
            if j == 0:
                ax.bar(xticks, this_sample_num, label=sample, width=0.6, edgecolor='w', color=SAMPLECOLOR[sample])
            else:
                bottom_sum = [sum(x) for x in zip(*bottoms[0: j])]
                ax.bar(xticks, this_sample_num, bottom=bottom_sum, label=sample, width=0.6, color=SAMPLECOLOR[sample])

        ax.set_xticks(xticks)
        ax.set_xticklabels(xlabels, fontsize=12, rotation=60)

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

        ax.set_yticks(np.linspace(0, 20, 5))
        ax.set_yticklabels([int(val) for val in np.linspace(0, 20, 5)], fontsize=12)

        if fig_idx == 0:
            ax.legend()
            ax.set_ylabel('Number', fontsize=13)
        fig_idx += 1

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/samples_somatic_mei_subtypes.pdf')
    plt.show()

if __name__ == '__main__':

    plot_ins_info('/data/home/jdlin/Prostate')