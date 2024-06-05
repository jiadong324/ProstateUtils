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


def dup_of_sample(workdir):

    dup_size = []

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


            if 'SVTYPE' in info_dict:
                svtype = info_dict['SVTYPE']
            else:
                svtype = entries[2].split('-')[-2]

            if svtype == 'DUP':
                dup_size.append([svlen, sample])

    print(len(dup_size))
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_ins_size = pd.DataFrame(dup_size, columns=['svlen', 'sample'])
    sns.histplot(data=df_ins_size, x='svlen', log_scale=True, hue='sample', alpha=1, bins=50,
                 hue_order=SAMPLES, palette=[SAMPLECOLOR[ele] for ele in SAMPLES], ax=ax)

    ax.set_yticks(np.linspace(0, 30, 4))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 30, 4)], fontsize=12)
    ax.set_ylabel('# of DUP', fontsize=13)

    ax.set_xlabel('SV length (bp)', fontsize=13)
    # ax.set_xticks([50, 100, 300, 1000, 6000, 10000, 100000])
    # ax.set_xticklabels(['50', r'$10^2$', '300', r'$10^3$', '6,000', r'$10^4$', '100,000'], fontsize=12, rotation=60)

    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=0)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/each_samples_dup_size.pdf')
    plt.show()

def samples_merged_dup(workdir):
    dup_size = []
    for line in open(f'{workdir}/sniffles_somatic/samples.somatic.merged.vcf', 'r'):
        if '#' in line:
            continue
        entries = line.strip().split('\t')
        info_dict = {}
        for token in entries[7].split(";"):
            if '=' in token:
                info_dict[token.split('=')[0]] = token.split('=')[1]

        svlen = abs(int(info_dict['SVLEN']))

        if 'SVTYPE' in info_dict:
            svtype = info_dict['SVTYPE']
        else:
            svtype = entries[2].split('-')[-2]
        on_gene = 'No'
        if svtype == 'DUP':
            if 'GENE' in info_dict:
                on_gene = 'Yes'
            dup_size.append([svlen, on_gene])

    print(len(dup_size))
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_dup_size = pd.DataFrame(dup_size, columns=['svlen', 'on_gene'])
    sns.histplot(data=df_dup_size, x='svlen', log_scale=True, hue='on_gene', hue_order=['Yes', 'No'], palette=['#fe9929', '#41b6c4'], alpha=1, bins=50, ax=ax)

    ax.set_yticks(np.linspace(0, 80, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 80, 5)], fontsize=12)
    ax.set_ylabel('# of DUP', fontsize=13)

    ax.set_xlabel('SV length (bp)', fontsize=13)

    ax.spines['top'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.tick_params(axis='both', which='both', length=0)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/all_samples_merged_dup_size.pdf')
    plt.show()

if __name__ == '__main__':
    # dup_of_sample('/data/home/jdlin/Prostate')
    samples_merged_dup('/data/home/jdlin/Prostate')