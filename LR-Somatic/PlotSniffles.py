#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/1

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
from matplotlib_venn import venn2

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



def somatic_brpk_counts(rev_v, aligner, min_supp):
    # svtypes = ['INS', 'DEL', 'DUP', 'INV', 'Inter-TRA', 'Intra-TRA']
    svtypes = ['INS', 'DEL', 'DUP', 'INV']

    svtypes_num = {svtype: [0 for i in range(len(SAMPLES))] for svtype in svtypes}
    # large_svtypes = {svtype: [0 for i in range(len(SAMPLES))] for svtype in ['INS', 'DEL', 'DUP', 'INV']}
    # sv_on_genes = {svtype: [0 for i in range(len(SAMPLES))] for svtype in svtypes}

    sv_num = {sample: [0 for i in range(len(VALID_CHROMS))] for sample in SAMPLES}
    sv_size = []
    insdel_size = {sample: [] for sample in SAMPLES}
    for i, sample in enumerate(SAMPLES):
        vcf_file = f'{REMOTEWORKDIR}/{sample}/sniffles/{rev_v}/{aligner}/{sample}.somatic.vcf'
        for line in open(vcf_file, 'r'):
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_dict = {}
            for token in entries[7].split(";"):
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            svtype = info_dict['SVTYPE']
            chr1, chr2 = entries[0], info_dict['CHR2']
            if svtype == 'BND':
                continue
            if int(info_dict['SUPPORT']) < min_supp:
                continue
            if abs(int(info_dict['SVLEN'])) > 100000:
                continue

            sv_size.append((sample, abs(int(info_dict['SVLEN'])), svtype))
            if svtype in ['INS', 'DEL']:
                insdel_size[sample].append((abs(int(info_dict['SVLEN'])), svtype))

            svtypes_num[svtype][i] += 1
            chrom_idx = VALID_CHROMS.index(chr1)
            sv_num[sample][chrom_idx] += 1

    xticks = np.arange(len(VALID_CHROMS))
    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    for sample, counts in sv_num.items():
        ax.plot(xticks, counts, color=SAMPLECOLOR[sample], lw=2, marker='o', label=sample)

    ax.legend()
    ax.set_xticks(xticks)
    ax.set_xticklabels(VALID_CHROMS, rotation=60, fontsize=13)
    ax.set_yticks(np.linspace(0, 400, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 400, 5)], fontsize=12)
    ax.set_ylabel('# of somatic SVs', fontsize=13)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/sniffles_somatic_bychrom.pdf')

    fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    xticks = np.arange((len(SAMPLES)))
    for j, svtype in enumerate(svtypes):
        this_data = svtypes_num[svtype]
        if j == 0:
            ax.bar(xticks, this_data, label=svtype, color=SVTYPECOLORS[svtype], width=0.8, edgecolor='w')
        else:
            bottom = []
            for m in range(0, j):
                bottom.append(svtypes_num[svtypes[m]])
            bottom_sum = [sum(x) for x in zip(*bottom)]
            ax.bar(xticks, this_data, bottom=bottom_sum, label=svtype, color=SVTYPECOLORS[svtype], width=0.8, edgecolor='w')

    ax.set_xticks(xticks)
    ax.set_xticklabels(SAMPLES, fontsize=13)
    ax.set_yticks(np.linspace(0, 2000, 5))
    ax.set_yticklabels([int(val) for val in np.linspace(0, 2000, 5)], fontsize=12)
    ax.legend()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylabel('# of somatic SVs', fontsize=13)
    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/sniffles_somatic_svtype_num.pdf')

    fig, ax = plt.subplots(1, 1, figsize=(7, 4))
    df_sv_size = pd.DataFrame(sv_size, columns=['sample', 'size', 'svtype'])
    sample_order = ['A', 'B', 'D', 'E', 'C']
    sample_color = [SAMPLECOLOR[ele] for ele in sample_order]
    legends = [Patch(label=sample, color=SAMPLECOLOR[sample]) for sample in SAMPLES]
    sns.histplot(data=df_sv_size, x='size', hue='sample', hue_order=sample_order, bins=50,
                 palette=sample_color, alpha=0.8, edgecolor='w', ax=ax, log_scale=True)

    ax.legend(handles=legends)
    ax.set_ylabel('# of somatic SVs', fontsize=13)
    ax.set_yticks(np.linspace(0, 400, 5))
    ax.set_yticklabels([int(ele) for ele in np.linspace(0, 400, 5)], fontsize=12)
    ax.set_xlabel('')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/sniffles_somatic_svsize.pdf')

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))
    df_insdel_C = pd.DataFrame(insdel_size['C'], columns=['size', 'svtype'])
    sns.histplot(data=df_insdel_C, x='size', hue='svtype', hue_order=['INS', 'DEL'], log_scale=True, bins=50,
                 palette=[SVTYPECOLORS['INS'], SVTYPECOLORS['DEL']], alpha=0.8, ax=ax)
    ax.set_ylabel('SV counts', fontsize=13)
    ax.set_yticks(np.linspace(0, 200, 5))
    ax.set_yticklabels([int(ele) for ele in np.linspace(0, 200, 5)], fontsize=12)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel('')
    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/t2t_based/Figure1/C_somatic_insdel_size.pdf')
    plt.show()

def somatic_brpk_sat(rev_v, aligner, min_supp):

    brpk_censat_num = {'Inside-Hsat':[0 for i in range(len(SAMPLES))],
                       'Outside-Hsat': [0 for i in range(len(SAMPLES))]}

    sample_brpk_censat = {sample: {} for sample in SAMPLES}
    sample_brkp_total = {sample: 0 for sample in SAMPLES}

    all_brpk = 0
    for i, sample in enumerate(SAMPLES):
        sample_brpk = 0
        vcf_file = f'{REMOTEWORKDIR}/{sample}/sniffles/{rev_v}/{aligner}/{sample}.somatic.vcf'
        for line in open(vcf_file, 'r'):
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_dict = {}
            for token in entries[7].split(";"):
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            svtype, sr = info_dict['SVTYPE'], int(info_dict['SUPPORT'])
            if svtype == 'BND':
                continue
            svlen = abs(int(info_dict['SVLEN']))
            if sr < min_supp and svlen >= 100000:
                continue

            chr1_censat = 'Outside-Hsat' if info_dict['SAT_POS1'] == 'None_censat' else 'Inside-Hsat'

            all_brpk += 1
            sample_brpk += 1
            if chr1_censat in brpk_censat_num:
                brpk_censat_num[chr1_censat][i] += 1
            else:
                brpk_censat_num[chr1_censat][i] = 1


            if info_dict['SAT_POS1'] in sample_brpk_censat[sample]:
                sample_brpk_censat[sample][info_dict['SAT_POS1']] += 1
            else:
                sample_brpk_censat[sample][info_dict['SAT_POS1']] = 1

            if svtype != 'INS':
                all_brpk += 1
                sample_brpk += 1
                chr2_censat = 'Outside-Hsat' if info_dict['SAT_POS2'] == 'None_censat' else 'Inside-Hsat'

                if chr2_censat in brpk_censat_num:
                    brpk_censat_num[chr2_censat][i] += 1
                else:
                    brpk_censat_num[chr2_censat][i] = 1

                if info_dict['SAT_POS2'] in sample_brpk_censat[sample]:
                    sample_brpk_censat[sample][info_dict['SAT_POS2']] += 1
                else:
                    sample_brpk_censat[sample][info_dict['SAT_POS2']] = 1


            # if svtype == 'BND':
            #     pos1_mindist_to_centro = get_brpk_to_centro(centro_pos[chr1], pos1)
            #     pos2_mindist_to_centro = get_brpk_to_centro(centro_pos[chr2], pos2)
            #     if chr1 == chr2:
            #         mindist_to_centro.append(min(pos1_mindist_to_centro, pos2_mindist_to_centro))
            #     else:
            #         mindist_to_centro.append(pos1_mindist_to_centro)
            #         mindist_to_centro.append(pos2_mindist_to_centro)
            # else:
            #     if svtype == 'INS':
            #         pos1_mindist_to_centro = get_brpk_to_centro(centro_pos[chr1], pos1)
            #         mindist_to_centro.append(pos1_mindist_to_centro)
            #     else:
            #         pos1_mindist_to_centro, pos2_mindist_to_centro = get_brpk_to_centro(centro_pos[chr1], pos1), get_brpk_to_centro(centro_pos[chr1], pos2)
            #         mindist_to_centro.append(min(pos1_mindist_to_centro, pos2_mindist_to_centro))
        sample_brkp_total[sample] = sample_brpk

    # sorted_brpk_censat = sorted(brpk_censat_num.items(), key=lambda x: x[1], reverse=True)[0: 10]
    # fig, ax = plt.subplots(1, 1, figsize=(5, 4))
    # xticks = np.arange(len(SAMPLES))
    # total = [i + j for i, j in zip(brpk_censat_num['Inside-Hsat'], brpk_censat_num['Outside-Hsat'])]
    # inside_hsat_pcrt = [i / j for i, j in zip(brpk_censat_num['Inside-Hsat'], total)]
    # outside_hsat_pcrt = [i / j for i, j in zip(brpk_censat_num['Outside-Hsat'], total)]
    #
    # ax.bar(xticks, inside_hsat_pcrt, width=0.7, edgecolor='w', color='#d95f02', label='Inside-HSat')
    # ax.bar(xticks, outside_hsat_pcrt, width=0.7, edgecolor='w', bottom=inside_hsat_pcrt, color='#1b9e77', label='Outside-HSat')
    # ax.set_xticks(xticks)
    # ax.set_xticklabels(SAMPLES, fontsize=13)
    #
    # ax.legend()
    # ax.set_ylim(0, 1)
    # ax.set_yticks(np.linspace(0, 1, 6))
    # ax.set_yticklabels([int(val * 100) for val in np.linspace(0, 1, 6)], fontsize=12)
    # ax.set_ylabel('Percent of somatic breakpoints', fontsize=13)
    #
    # ax.spines['left'].set_visible(False)
    # ax.spines['right'].set_visible(False)
    # fig.tight_layout()

    fig, axes = plt.subplots(1, 5, figsize=(16, 4))
    yticks = np.arange(5)
    for sample, censat_annot in sample_brpk_censat.items():
        sorted_censat_annot = sorted(censat_annot.items(), key=lambda x:x[1], reverse=True)[0: 5]
        fig_idx = SAMPLES.index(sample)
        axes[fig_idx].set_title(f'{sample}', fontsize=13)
        axes[fig_idx].barh(yticks, [ele[1] / sample_brkp_total[sample] for ele in sorted_censat_annot], height=0.5)
        axes[fig_idx].set_yticks(yticks)
        axes[fig_idx].set_yticklabels([ele[0] for ele in sorted_censat_annot])
        axes[fig_idx].spines['top'].set_visible(False)
        axes[fig_idx].spines['right'].set_visible(False)
        axes[fig_idx].invert_yaxis()

    fig.tight_layout()
    plt.show()


def somatic_insdel_size(ref_v, min_supp, aligner):
    # svtypes = ['INS', 'DEL', 'DUP', 'INV']
    svtypes = ['INS', 'DEL']
    svsize = {ele: [] for ele in svtypes}

    for i, sample in enumerate(SAMPLES):
        vcf_file = f'{REMOTEWORKDIR}/{sample}/sniffles/{ref_v}/{aligner}/{sample}.somatic.vcf'
        for line in open(vcf_file, 'r'):
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_dict = {}
            for token in entries[7].split(";"):
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]
            svtype = info_dict['SVTYPE']

            if svtype in svtypes:
                sr, svlen = int(info_dict['SUPPORT']), abs(int(info_dict['SVLEN']))
                if sr >= min_supp and svlen < 100000:
                    svsize[svtype].append((svlen, sample))

    sample_order = ['A', 'B', 'D', 'E', 'C']
    sample_color = [SAMPLECOLOR[ele] for ele in sample_order]
    legends = [Patch(label=sample, color=SAMPLECOLOR[sample]) for sample in SAMPLES]

    for svtype, sizes in svsize.items():
        df_size = pd.DataFrame(sizes, columns=['size', 'sample'])
        fig, ax = plt.subplots(1, 1, figsize=(6, 4))
        sns.histplot(data=df_size, hue='sample', x='size', log_scale=True, alpha=0.6, bins=50,
                     hue_order=sample_order, palette=sample_color, ax=ax)

        # medians = df_size.groupby(['svtype'])['size'].median()

        ax.set_ylabel(f'# of {svtype}', fontsize=12)
        ax.set_xlabel('')

        # for xtick in ax.get_xticks():
        #     ax.text(xtick, round(medians[xtick], 2), round(medians[xtick], 2),
        #                   horizontalalignment='center', size='x-small', color='k', weight='semibold')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)

        ax.set_xticks([50, 100, 300, 1000, 10000, 100000])
        ax.set_xticklabels(['50', '100', '300', '1,000', '10,000', '100,000'], fontsize=13, rotation=60)
        ax.legend(handles=legends)
        fig.tight_layout()
        # fig.savefig(f'{REMOTEFIG}/{sample}_aligners_merged_somatic_svsize.pdf')

    plt.show()

def somatic_insdel_rep(ref_v, aligner, min_supp):

    # rep_class = {ele: [0 for i in range(len(SAMPLES))] for ele in GENOMICREGIONS}
    rep_class = {sample: {} for sample in SAMPLES}
    sample_svs = {sample: 0 for sample in SAMPLES}
    for i, sample in enumerate(SAMPLES):
        sample_brpk = 0
        vcf_file = f'{REMOTEWORKDIR}/{sample}/sniffles/{ref_v}/{aligner}/{sample}.somatic.vcf'
        for line in open(vcf_file, 'r'):
            if '#' in line:
                continue
            entries = line.strip().split('\t')
            info_dict = {}
            for token in entries[7].split(";"):
                if '=' in token:
                    info_dict[token.split('=')[0]] = token.split('=')[1]

            svtype, sr = info_dict['SVTYPE'], int(info_dict['SUPPORT'])
            if svtype == 'BND':
                continue
            svlen = abs(int(info_dict['SVLEN']))
            if sr < min_supp and svlen >= 100000:
                continue
            sample_brpk += 1
            chr1_class, chr1_subclass = info_dict['REP_POS1'].split(',')

            if chr1_class in rep_class[sample]:
                rep_class[sample][chr1_class] += 1
            else:
                rep_class[sample][chr1_class] = 1

            if svtype != 'INS':
                sample_brpk += 1
                chr2_class, chr2_subclass = info_dict['REP_POS2'].split(',')

                if chr2_class in rep_class[sample]:
                    rep_class[sample][chr2_class] += 1
                else:
                    rep_class[sample][chr2_class] = 1
        sample_svs[sample] = sample_brpk

    fig, axes = plt.subplots(1, 5, figsize=(16, 4))
    yticks = np.arange(5)
    for sample, reps in rep_class.items():
        fig_idx = SAMPLES.index(sample)
        sorted_reps = sorted(reps.items(), key=lambda x: x[1], reverse=True)[0:5]
        axes[fig_idx].barh(yticks, [ele[1] / sample_svs[sample] for ele in sorted_reps], height=0.5)
        axes[fig_idx].set_title(f'{sample}', fontsize=13)
        axes[fig_idx].set_yticks(yticks)
        axes[fig_idx].set_yticklabels([ele[0] for ele in sorted_reps])
        axes[fig_idx].spines['top'].set_visible(False)
        axes[fig_idx].spines['right'].set_visible(False)
        axes[fig_idx].invert_yaxis()

    fig.tight_layout()
    plt.show()


def make_autopct(values):
    def my_autopct(pct):
        total = sum(values)
        val = int(round(pct*total/100.0))
        return '{p:.2f}%  ({v:d})'.format(p=pct,v=val)
    return my_autopct

def get_brpk_to_centro(centro_pos, pos):
    dists = []
    for (start, end) in centro_pos:
        pos_dist1, pos_dist2 = abs(pos - start), abs(pos - end)
        dists.extend([pos_dist1, pos_dist2])

    return sorted(dists)[0]



def main():
    min_supp = 6

    somatic_brpk_counts('T2T', 'minimap2', min_supp)
    # somatic_insdel_rep('T2T', 'minimap2', min_supp)
    # somatic_brpk_sat('T2T', 'minimap2', min_supp)

    # somatic_insdel_size('T2T', min_supp, 'minimap2')



if __name__ == '__main__':
    main()

'''

def plot_sv_bp_density():
    window_size = 10000000

    chrom_size_dict = {}
    chrom_binned_trees = {}
    chrom_binned_list = []
    most_bins = 0
    with open(CHROMSIZE, 'r') as f:
        for line in f:
            entries = line.strip().split('\t')
            chrom, size = entries[0], int(entries[1])
            if chrom in VALID_CHROMS:
                chrom_size_dict[chrom] = size

                chrom_binned_trees[chrom] = IntervalTree()
                bins = int(size / window_size)
                start = 0
                for ith in range(bins):
                    end = start + window_size
                    chrom_binned_trees[chrom][start: end] = (ith, start)
                    start = end

                if start < size:
                    chrom_binned_trees[chrom][start: size] = (bins, start)

                if chrom == 'chr1':
                    most_bins = bins + 1
                chrom_binned_list.append([0 for i in range(most_bins)])

    columns = ['chrom', 'start', 'end', 'id', 'svtype', 'svlen']
    for sample in SAMPLES:
        all_sv_bp = pd.read_csv(f'{SOMATICDIR}/merged_aligner/{sample}.minimap2.ngmlr.merged.bed', sep='\t',
                                names=columns)

        for idx, row in all_sv_bp.iterrows():
            chrom, start, end = row['chrom'], int(row['start']), int(row['end'])

            overlaps = list(chrom_binned_trees[chrom].overlap(start, end))
            contig_idx = VALID_CHROMS.index(chrom)
            for overlap in overlaps:
                bin_num = overlap.data[0]
                chrom_binned_list[contig_idx][bin_num] += 1

    fig, ax = plt.subplots(1, 1, figsize=(6, 4))

    sns.heatmap(chrom_binned_list, ax=ax, cmap="YlGnBu")
    ax.set_yticks([r + 0.5 for r in np.arange(len(VALID_CHROMS))])
    ax.set_yticklabels(VALID_CHROMS, rotation=360)

    ax.set_xticks([r + 0.5 for r in np.arange(25)])
    ax.set_xticklabels([int(val) for val in np.arange(25)])

    fig.tight_layout()
    plt.show()
    fig.savefig(f'{SOMATICFIGURES}/breakpoint_density.pdf')

def aligner_concordance(workdir):

    fig, axes = plt.subplots(1, len(SAMPLES), figsize=(15, 5))

    for i, sample in enumerate(SAMPLES):
        ax = axes[i]
        ax.set_title(f'Sample: {sample}', fontsize=13)
        aligner_merged_vcf = f'{workdir}/{sample}/sniffles/GRCh38/{sample}.somatic.aligners.jasmine.vcf'
        supp_vec, merged_total = get_survivor_suppvec(aligner_merged_vcf)
        venn2(subsets={'10': supp_vec['10'], '01': supp_vec['01'], '11': supp_vec['11']}, set_labels=['ngmlr', 'minimap2'], ax=ax)

    fig.tight_layout()
    fig.savefig(f'{REMOTEFIG}/sample_aligner_overlap_somatic_svs.pdf')
    plt.show()


def samples_sv_upsetplot(workdir):

    merged_vcf = f'{workdir}/sniffles_somatic/samples.somatic.merged.vcf'
    suppvec_dict, merged_total = get_survivor_suppvec(merged_vcf)

    upset_mem = []
    data = []
    for suppvec, count in suppvec_dict.items():
        this_mem = []
        for i, val in enumerate(suppvec):
            if val == '1':
                this_mem.append(SAMPLES[i])
        upset_mem.append(this_mem)
        data.append(count)
    print(data)
    upset_data = upsetplot.from_memberships(upset_mem, data=data)
    upsetplot.UpSet(upset_data).plot()
    plt.savefig(f'{REMOTEFIG}/samples_somatic_sv_upsetplot.pdf')
    plt.show()

'''
