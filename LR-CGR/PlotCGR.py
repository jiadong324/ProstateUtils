#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/6/30

'''

import matplotlib.pylab as plt
import seaborn as sns
import numpy as np
from matplotlib.patches import Patch


from Helpers.Constants import *

plt.rcParams["pdf.fonttype"] = 42
plt.rcParams["xtick.major.width"] = 2
plt.rcParams["xtick.minor.width"] = 2
plt.rcParams["xtick.labelsize"] = 11

plt.rcParams["ytick.major.width"] = 2
plt.rcParams["ytick.minor.width"] = 2
plt.rcParams["ytick.labelsize"] = 11

plt.rcParams["axes.linewidth"] = 2


def cgr_reads_brpk_counts(aligner):
    reads_brpk = []
    for j, sample in enumerate(SAMPLES):
        for line in open(f'{REMOTEWORKDIR}/{sample}/somatic_trans/T2T/{aligner}/all_sa/somatic_trans_junctions_sreads.bed', 'r'):
            entries = line.strip().split('\t')
            reads_brpk.append((sample, int(entries[1])))

def cgr_brpk_counts(aligner, min_sr):

    brpk_nums = {'Intra': [0 for i in range(len(SAMPLES))], 'Inter': [0 for i in range(len(SAMPLES))]}
    link_size_dict = {}
    for j, sample in enumerate(SAMPLES):
        for line in open(f'{REMOTEWORKDIR}/{sample}/somatic_trans/T2T/{aligner}/all_sa/somatic_trans_junctions_links.tsv', 'r'):
            entries = line.strip().split('\t')
            link_id, supp, brpks_info, brpks_af, brpk_censat = entries[0], int(entries[1]), entries[2].split(';'), entries[4].split(';'), entries[6].split(';')
            link_size = len(link_id.split('-'))
            if supp < min_sr:
                continue

            if link_size in link_size_dict:
                link_size_dict[link_size] += 1
            else:
                link_size_dict[link_size] = 1

            chrom_set = set()
            for i in range(len(brpks_info)):
                chrom, pos = brpks_info[i].split(':')[0], int(brpks_info[i].split(':')[1])
                chrom_set.add(chrom)

            if len(chrom_set) == 1:
                brpk_nums['Intra'][j] += 1
            else:
                brpk_nums['Inter'][j] += 1

    # legends = [Patch(color='#1b9e77', label='Intra'), Patch(color='#d95f02', label='Inter')]
    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    xticks = np.arange(len(SAMPLES))

    ax.bar(xticks, brpk_nums['Intra'], color=SVTYPECOLORS['Intra-TRA'], edgecolor='w', label='Intra-chromosome')
    ax.bar(xticks, brpk_nums['Inter'], bottom=brpk_nums['Intra'], color=SVTYPECOLORS['Inter-TRA'], edgecolor='w', label='Inter-chromosome')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    ax.set_yticks(np.linspace(0, 300, 4))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 300, 4)])
    ax.set_ylabel('#TRA breakpoints', fontsize=13)

    ax.legend(fontsize=13, loc='lower right')

    ax.set_xticks(xticks)
    ax.set_xticklabels(SAMPLES, fontsize=13)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/tra_breakpoints_counts.pdf')

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    sorted_link_size = sorted(link_size_dict.items(), key=lambda x:x[1], reverse=True)
    values = []
    labels = []

    larg = 0
    for (key, val) in sorted_link_size:
        if key < 4:
            values.append(val)
            labels.append(key)
        else:
            larg += val
    values.append(larg)
    labels.append('>=4')
    print(values)

    ax.pie(values, labels=labels, autopct=make_autopct(values), colors=['#dba987','#a3c1b1', '#adaccc'], explode=[0, 0, 0.2])
    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/spanned_tra_breakpoints.pdf')
    plt.show()


def cgr_brpk_rep(aligner, min_sr):

    brpk_by_chrom = {chrom: {} for chrom in VALID_CHROMS}
    brpk_num = {chrom: 0 for chrom in VALID_CHROMS}
    brpk_censat_num = {}
    all_brpk = 0
    brpk_pos = {chrom: [] for chrom in VALID_CHROMS}
    for sample in SAMPLES:

        for line in open(f'{REMOTEWORKDIR}/{sample}/somatic_trans/T2T/{aligner}/all_sa/somatic_trans_junctions_links.tsv', 'r'):
            entries = line.strip().split('\t')
            link_id, supp, brpks_info, brpks_af, brpk_censat = entries[0], int(entries[1]), entries[2].split(';'), entries[5].split(';'), entries[7].split(';')
            brpk_mapq = entries[3].split(';')
            if supp < min_sr:
                continue
            for i in range(len(brpks_info)):
                chrom, pos, mapq = brpks_info[i].split(':')[0], int(brpks_info[i].split(':')[1]), int(brpk_mapq[i])
                censat = brpk_censat[i]
                brpk_num[chrom] += 1

                all_brpk += 1
                if censat in brpk_by_chrom[chrom]:
                    brpk_by_chrom[chrom][censat] += 1
                else:
                    brpk_by_chrom[chrom][censat] = 1
                brpk_pos[chrom].append(pos / 100000)

                if censat in brpk_censat_num:
                    brpk_censat_num[censat] += 1
                else:
                    brpk_censat_num[censat] = 1

    sorted_brpk_censat = sorted(brpk_censat_num.items(), key=lambda x:x[1], reverse=True)[0: 10]

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    pie_values = [sorted_brpk_censat[0][1], all_brpk - sorted_brpk_censat[0][1]]
    ax.pie(pie_values, autopct=make_autopct(pie_values), colors=['#cbce6b', '#8dbfc9'], explode=[0, 0.1], startangle=120)
    legends = [Patch(label='Outside Hsat', color='#cbce6b'), Patch(label='Inside Hsat', color='#8dbfc9')]
    ax.legend(handles=legends, loc='upper left')
    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/tra_breakpoints_hsat.pdf')

    # fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    # censat_brpk = sorted(brpk_censat_num.items(), key=lambda x:x[1], reverse=True)[1:5]
    #
    # xticks = np.arange(len(censat_brpk))
    #
    # ax.bar(xticks, [ele[1]/ (all_brpk - sorted_brpk_censat[0][1]) for ele in censat_brpk])
    # ax.set_xticks(xticks)
    # ax.set_xticklabels([ele[0] for ele in censat_brpk], rotation=60)
    #
    # ax.spines['left'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    #
    # fig.tight_layout()

    # for chrom, brpk_censat_count in brpk_by_chrom.items():
    #     if chrom == 'chr9':
    #         fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    #         ax.set_title(chrom, fontsize=13)
    #         sorted_brpk_censat_count = sorted(brpk_censat_count.items(), key=lambda x:x[1], reverse=True)[0:3]
    #         brpk_pcrt = [ele[1] / brpk_num[chrom] for ele in sorted_brpk_censat_count]
    #
    #         xticks = np.arange(len(sorted_brpk_censat_count))
    #
    #         bars = ax.bar(xticks, brpk_pcrt, width=0.5)
    #         # bars[0].set_color('#e7298a')
    #
    #         ax.spines['left'].set_visible(False)
    #         ax.spines['right'].set_visible(False)
    #
    #         ax.set_xticks(xticks)
    #         ax.set_xticklabels([ele[0] for ele in sorted_brpk_censat_count], rotation=60)
    #         fig.tight_layout()
    #
    #         fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    #         yvals = [1 for i in range(len(brpk_pos[chrom]))]
    #         ax.axvline(44951775 / 100000, color='r')
    #         ax.axvline(47582595 / 100000, color='r')
    #         sns.histplot(brpk_pos[chrom], ax=ax)
    #         ax.spines['left'].set_visible(False)
    #         ax.spines['right'].set_visible(False)
    #         fig.tight_layout()

    plt.show()

def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def main():

    min_supp = 6
    aligner = 'minimap2'

    cgr_brpk_counts(aligner, min_supp)
    # cgr_brpk_rep(aligner, min_supp)

if __name__ == '__main__':
    main()