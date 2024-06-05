#!/usr/bin/env python

# encoding: utf-8

'''

@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/12/29
'''

import math

import matplotlib.pylab as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
import pysam
import seaborn as sns
from matplotlib.patches import Patch


from Helpers.Constants import *
from Helpers.Functions import *


plt.rcParams["font.family"] = "Arial"
plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2



def dv_phase_stats(workdir):

    n50_length = []
    total_variants = []
    het_snvs = []
    phased_pcrt = []
    chroms_n50 = []
    chroms = set()

    for sample in SAMPLES:
        df_phase = pd.read_csv(f'{workdir}/{sample}/phased_snps/dv/{sample}-T.phased.stats', sep='\t', header=[0])
        longest_n50 = [0, '']
        for idx, row in df_phase.iterrows():
            variants, het_snv, phased, unphased, n50 = int(row['variants']), int(row['heterozygous_snvs']), int(row['phased_snvs']), int(row['unphased']), int(row['block_n50']) / 1000
            if row['chromosome'] == 'ALL':
                n50_length.append(n50)
                total_variants.append(variants)
                het_snvs.append([sample, het_snv / 1000, 'het_total'])
                het_snvs.append([sample, phased / 1000, 'het_phased'])
                phased_pcrt.append(phased / het_snv * 100)

            elif row['chromosome'] != 'chrX':
                chroms_n50.append([n50, sample, row['chromosome']])
                if n50 > longest_n50[0]:
                    longest_n50 = [n50, row['chromosome']]
                chroms.add(row['chromosome'])
        print(f'{sample}, longest N50 block: {longest_n50[1]}, {longest_n50[0]}')

    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    xticks = np.arange(len(chroms))

    df_chroms_n50 = pd.DataFrame(chroms_n50, columns=['n50', 'sample', 'chrom'])
    sns.lineplot(data=df_chroms_n50, x='chrom', hue='sample', y='n50', marker='o', lw=2)

    ax.set_ylim(00, 400)
    ax.set_yticks(np.linspace(00, 400, 5))
    ax.set_yticklabels([0, 100, 200, 300, 400], fontsize=12)
    ax.set_ylabel('N50 phase block (Kbp)', fontsize=13)

    ax.set_xticks(xticks)
    ax.set_xticklabels(VALID_CHROMS[0: 22], fontsize=13, rotation=60)
    ax.set_xlabel('')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    # fig.savefig(f'{REMOTEFIG}/samples_dv_n50_block.pdf')


    df_het_snps = pd.DataFrame(het_snvs, columns=['sample', 'num', 'status'])
    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    sns.barplot(data=df_het_snps, x='sample', y='num', hue='status', ax=ax)

    ax.set_xticks(np.arange(len(SAMPLES)))
    ax.set_xticklabels(SAMPLES, fontsize=13)
    ax.set_xlabel('')
    ax.legend(loc='lower left')

    ax.set_ylabel('Count (x1000)', fontsize=13)
    ax.set_yticks(np.linspace(0, 2000, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 2000, 5)], fontsize=12)

    ax.spines['top'].set_visible(False)

    ax_twin = ax.twinx()
    ax_twin.plot(np.arange(len(SAMPLES)), phased_pcrt, lw=1.5, color='r', marker='o')
    ax_twin.set_ylabel('Percentage of phased', fontsize=13)
    ax_twin.set_yticks(np.linspace(0, 100, 5))
    ax_twin.set_yticklabels([int(ele) for ele in np.linspace(0, 100, 5)], fontsize=12)

    ax_twin.spines['top'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/samples_dv_het_snps.pdf')
    plt.show()

# def phased_somatic_snvs(workdir):



def main():

    dv_phase_stats('/data/home/jdlin/Prostate')
    # varscan_phase_stats('/data/home/jdlin/Prostate')

if __name__ == '__main__':
    main()

