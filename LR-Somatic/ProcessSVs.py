#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/3/1

'''

import vcf
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import os
import sys

sys.path.append("//")

from Helpers.Constants import *
from Helpers.Functions import *



def separate_sniffles_somatic(workdir, ref_v, aligner):
    '''
    Separate somatic SVs from tumor-normal paired calling
    :param workdir:
    :param ref_v:
    :return:
    '''

    cosmic_region, region_gene_map = parse_cosmic_gene()


    # simple_reps = pysam.Tabixfile(REMOTESIM, 'r')
    # rmsk = pysam.Tabixfile(REMOTERMSK, 'r')
    # sds = pysam.Tabixfile(REMOTESD, 'r')

    rmsk = pysam.Tabixfile('/data/home/jdlin/reference/T2T/hs1.repeatMasker.sorted.bed.gz', 'r')
    censat = pysam.Tabixfile('/data/home/jdlin/reference/T2T/chm13v2.0_censat_v2.0.bed.gz', 'r')

    somatic_sv_counts = open(f'{workdir}/{ref_v}_{aligner}_somatic_nums.tsv', 'w')

    for sample in SAMPLES:
        somatic_bnd_writer = open(f'{workdir}/{sample}/sniffles/{ref_v}/{aligner}/{sample}.somatic.bnd.tsv', 'w')

        vcf_input = open(f'{workdir}/{sample}/sniffles/{ref_v}/{aligner}/{sample}.tumor-normal.vcf', 'r')
        # somatic_exbnd_writer = open(f'{workdir}/{sample}/{caller}/{ref_v}/{sample}.somatic.exbnd.bed', 'w')
        somatic_vcf_writer = open(f'{workdir}/{sample}/sniffles/{ref_v}/{aligner}/{sample}.somatic.vcf', 'w')
        # somatic_ex_rep_vcf = open(f'{workdir}/{sample}/sniffles/{ref_v}/{aligner}/{sample}.somatic.exrep.vcf', 'w')

        sv_count_dict = {'INS': 0, 'DEL': 0, 'DUP': 0, 'BND': 0, 'INV': 0}
        all_somatic_num = 0

        for line in vcf_input:

            if '#' in line:
                print(line.strip(), file=somatic_vcf_writer)
                # print(line.strip(), file=somatic_ex_rep_vcf)
                continue

            entries = line.strip().split('\t')
            chrom, start, id = entries[0], int(entries[1]), entries[2]

            if chrom not in VALID_CHROMS:
                continue

            info_tokens = entries[7].split(';')
            if info_tokens[0] == 'IMPRECISE':
                continue

            info_dict = {}
            for token in info_tokens[1:]:
                info_dict[token.split('=')[0]] = token.split('=')[1]

            sup_vec, svtype, sr = info_dict['SUPP_VEC'], info_dict['SVTYPE'], int(info_dict['SUPPORT'])
            info_str1 = '\t'.join(entries[0:7])
            info_str2 = '\t'.join(entries[8:])

            if sr < 6:
                continue

            # pos1_region_label, pos1_rptype, pos1_pcrt = annotate_sv_region(chrom, start - 10, start + 10, 0, simple_reps, rmsk, sds)
            start_iv = cosmic_region[chrom].overlap(start - 50, start + 50)
            # if start_iv:

            pos1_repclass, pos1_repsubtype, pos1_pcrt = rmsk_annotation(chrom, start - 10, start + 10, 0, rmsk)
            pos1_cen, pos1_pcrt = censat_annotation(chrom, start - 10, start + 10, 0, censat)

            if sup_vec == '01':
                if svtype == 'BND':
                    chr2, pos2 = get_bnd_contigs(entries[4])
                    if chr2 in VALID_CHROMS:

                        # pos2_region_label, pos2_rptype, pos2_pcrt = annotate_sv_region(chr2, pos2 - 10, pos2 + 10, 0, simple_reps, rmsk, sds)
                        end_iv = cosmic_region[chrom].overlap(pos2 - 50, pos2 + 50)
                        pos2_repclass, pos2_repsubtype, pos2_pcrt = rmsk_annotation(chr2, pos2 - 10, pos2 + 10, 0, rmsk)
                        pos2_cen, pos2_pcrt = censat_annotation(chr2, pos2 - 10, pos2 + 10, 0, censat)

                        print(f'{chrom}\t{start}\t{pos1_repclass}\t{chr2}\t{pos2}\t{pos2_repclass}\t{sr}\t{aligner}', file=somatic_bnd_writer)

                        new_info_str = f'{entries[7]};CHR2={chr2};END={pos2};REP_POS1={pos1_repclass},{pos1_repsubtype};' \
                                       f'REP_POS2={pos2_repclass},{pos2_repsubtype};SAT_POS1={pos1_cen};SAT_POS2={pos2_cen}'

                        print(f'{info_str1}\t{new_info_str}\t{info_str2}', file=somatic_vcf_writer)
                        sv_count_dict['BND'] += 1

                        continue

                else:
                    sv_count_dict[svtype] += 1
                    all_somatic_num += 1
                    svlen = int(info_dict['SVLEN'])
                    end = start + svlen if svtype == 'INS' else int(info_dict['END'])
                    if svtype != 'INS':
                        # pos2_region_label, pos2_rptype, pos2_pcrt = annotate_sv_region(chrom, end - 10, end + 10, 0, simple_reps, rmsk, sds)
                        pos2_repclass, pos2_repsubtype, pos2_pcrt = rmsk_annotation(chrom, end - 10, end + 10, 0, rmsk)
                        pos2_cen, pos2_pcrt = censat_annotation(chrom, end - 10, end + 10, 0, censat)

                        new_info_str = f'{entries[7]};CHR2={chrom};REP_POS1={pos1_repclass},{pos1_repsubtype};' \
                                       f'REP_POS2={pos2_repclass},{pos2_repsubtype};SAT_POS1={pos1_cen};SAT_POS2={pos2_cen}'
                        print(f'{info_str1}\t{new_info_str}\t{info_str2}', file=somatic_vcf_writer)
                    else:
                        new_info_str = f'{entries[7]};CHR2={chrom};REP_POS1={pos1_repclass},{pos1_repsubtype};SAT_POS1={pos1_cen}'
                        print(f'{info_str1}\t{new_info_str}\t{info_str2}', file=somatic_vcf_writer)

        out_str = str(sv_count_dict['INS']) + '\t' + str(sv_count_dict['DEL']) + '\t' + str(sv_count_dict['DUP']) + '\t' + str(sv_count_dict['INV']) + '\t' + str(sv_count_dict['BND'])
        print(f'{sample}\t{aligner}\t{all_somatic_num}\t{out_str}', file=somatic_sv_counts)

def add_haptag(workdir):
    '''
    Add haplotag to TRA detected via the graph-based approach and sniffles
    :param workdir:
    :return:
    '''
    for sample in SAMPLES:
        phased_reads_dict = {}
        phased_reads = pd.read_csv(f'{workdir}/{sample}.phased.reads.txt', sep='\t', header=[0])

        for idx, row in phased_reads.iterrows():
            qname = row['#readname']
            if 'ccs' in qname:
                phased_reads_dict[qname] = row['haplotype']

        # haptagged_trans = open(f'{workdir}/{sample}_somatic_trans_junctions_links.haptag.tsv', 'w')
        # for line in open(f'{workdir}/{sample}_somatic_trans_junctions_links.tsv', 'r'):
        #     entries = line.strip().split('\t')
        #     haptag = ''
        #     for read in entries[8].split(';'):
        #         if read in phased_reads_dict:
        #             haptag += str(phased_reads_dict[read])
        #     if haptag == '':
        #         print(f'{line.strip()}\tUnphased',file=haptagged_trans)
        #     else:
        #         print(f'{line.strip()}\t{haptag}', file=haptagged_trans)

        haptag_svs = open(f'{workdir}/{sample}.somatic.haptag.vcf', 'w')
        bnd_bed = open(f'{workdir}/{sample}.somatic.BND.bed', 'w')
        for line in open(f'{workdir}/{sample}.somatic.vcf', 'r'):
            if '#' in line:
                print(line.strip(), file=haptag_svs)
                continue
            entries = line.strip().split('\t')
            info_tokens = entries[7].split(";")
            info_dict = {}

            for token in info_tokens:
                if "=" not in token:
                    continue
                info_dict[token.split("=")[0]] = token.split("=")[1].replace(">", "")
            haptag = ''
            for read in info_dict['RNAMES'].split(','):
                if read in phased_reads_dict:
                    haptag += str(phased_reads_dict[read])
            str1 = '\t'.join(entries[0:7])
            str2 = '\t'.join(entries[8:])
            new_info_str = entries[7] + f';HAP=-1' if haptag == '' else entries[7] + f';HAP={haptag}'

            if info_dict['SVTYPE'] == 'BND':
                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(entries[0], entries[1], int(entries[1]) + 1, info_dict['CHR2'], info_dict['END'],
                                                            int(info_dict['END']) + 1, info_dict['SUPPORT']), file=bnd_bed)
            elif int(info_dict['SVLEN']) >= 100000:
                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}'.format(entries[0], entries[1], int(entries[1]) + 1,
                                                             entries[0], info_dict['END'],
                                                             int(info_dict['END']) + 1, info_dict['SUPPORT']), file=bnd_bed)

            new_vcf_str = f'{str1}\t{new_info_str}\t{str2}'
            print(new_vcf_str, file=haptag_svs)




def annot_severus(workdir):

    rmsk = pysam.Tabixfile('/Users/jiadonglin/Data/genome/T2T/hs1.repeatMasker.sorted.bed.gz', 'r')
    censat = pysam.Tabixfile('/Users/jiadonglin/Data/genome/T2T/chm13v2.0_censat_v2.0.bed.gz', 'r')

    for sample in SAMPLES:
        annot_vcf = open(f'{workdir}/{sample}.severus_somatic.annot.vcf', 'w')
        bnd_bed = open(f'{workdir}/{sample}.severus_somatic.BND.bed', 'w')

        for line in open(f'{workdir}/{sample}.severus_somatic.vcf'):
            if '#' in line:
                print(line.strip(), file=annot_vcf)
                continue

            entries = line.strip().split('\t')
            chrom, start, id = entries[0], int(entries[1]), entries[2]

            if chrom not in VALID_CHROMS:
                continue

            info_tokens = entries[7].split(';')
            if info_tokens[0] == 'IMPRECISE':
                continue
            info_str1 = '\t'.join(entries[0:7])
            info_str2 = '\t'.join(entries[8:])

            info_dict = {}
            for token in info_tokens[1:]:
                info_dict[token.split('=')[0]] = token.split('=')[1]
            svtype = info_dict['SVTYPE']
            pos1_repclass, pos1_repsubtype, pos1_pcrt = rmsk_annotation(chrom, start - 10, start + 10, 0, rmsk)
            pos1_cen, pos1_pcrt = censat_annotation(chrom, start - 10, start + 10, 0, censat)


            if svtype == 'BND':
                chr2, pos2 = info_dict['CHR2'], int(info_dict['END'])
                if chr2 not in VALID_CHROMS:
                    continue
                pos2_repclass, pos2_repsubtype, pos2_pcrt = rmsk_annotation(chr2, pos2 - 10, pos2 + 10, 0, rmsk)
                pos2_cen, pos2_pcrt = censat_annotation(chr2, pos2 - 10, pos2 + 10, 0, censat)
                new_info_str = f'{entries[7]};CHR2={chr2};END={pos2};REP_POS1={pos1_repclass},{pos1_repsubtype};' \
                               f'REP_POS2={pos2_repclass},{pos2_repsubtype};SAT_POS1={pos1_cen};SAT_POS2={pos2_cen}'

                print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(entries[0], entries[1], int(entries[1]) + 1,
                                                             info_dict['CHR2'], info_dict['END'],
                                                             int(info_dict['END']) + 1), file=bnd_bed)
                print(f'{info_str1}\t{new_info_str}\t{info_str2}', file=annot_vcf)
            else:
                end = int(info_dict['END'])
                if int(info_dict['SVLEN']) >= 100000:
                    print('{0}\t{1}\t{2}\t{3}\t{4}\t{5}'.format(entries[0], entries[1], int(entries[1]) + 1,
                                                                 entries[0], info_dict['END'],
                                                                 int(info_dict['END']) + 1), file=bnd_bed)
                    continue

                if svtype != 'INS':
                    pos2_repclass, pos2_repsubtype, pos2_pcrt = rmsk_annotation(chrom, end - 10, end + 10, 0, rmsk)
                    pos2_cen, pos2_pcrt = censat_annotation(chrom, end - 10, end + 10, 0, censat)
                    new_info_str = f'{entries[7]};CHR2={chrom};REP_POS1={pos1_repclass},{pos1_repsubtype};' \
                                   f'REP_POS2={pos2_repclass},{pos2_repsubtype};SAT_POS1={pos1_cen};SAT_POS2={pos2_cen}'
                    print(f'{info_str1}\t{new_info_str}\t{info_str2}', file=annot_vcf)
                else:
                    new_info_str = f'{entries[7]};CHR2={chrom};REP_POS1={pos1_repclass},{pos1_repsubtype};SAT_POS1={pos1_cen}'
                    print(f'{info_str1}\t{new_info_str}\t{info_str2}', file=annot_vcf)

def censat_annotation(chrom, start, end, pcrt_thresh, censat_tabix):
    if start > end:
        start, end = end, start
    size = end - start + 1

    annotations = []
    for censat in censat_tabix.fetch(chrom, start, end):
        entries = censat.strip().split('\t')
        rp_start, rp_end, subtype = int(entries[1]), int(entries[2]), entries[3]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            if overlap_pcrt >= pcrt_thresh:
                annotations.append((subtype, overlap_pcrt))

    if len(annotations) == 0:
        return ('None_censat', 0)

    sorted_annotations = sorted(annotations, key=lambda x: x[1], reverse=True)
    return sorted_annotations[0]

def rmsk_annotation(chrom, start, end, pcrt_thresh, rmsk_tabix):

    if start > end:
        start, end = end, start
    size = end - start + 1
    annotations = []
    for rmsk in rmsk_tabix.fetch(chrom, start, end):
        entries = rmsk.strip().split('\t')
        rp_start, rp_end, subtype, repclass = int(entries[1]), int(entries[2]), entries[3], entries[4]
        overlap_size = get_overlaps(start, end, rp_start, rp_end)
        if overlap_size > 0:
            overlap_pcrt = min(overlap_size / size * 100, 100)
            if overlap_pcrt >= pcrt_thresh:
                if repclass == 'Simple_repeat':
                    motif = subtype[1: -2]
                    new_subtype = 'VNTR' if len(motif) >= 7 else 'STR'
                    annotations.append((new_subtype, subtype, overlap_pcrt))
                    # return ('Tandem Repeats', subtype, overlap_pcrt)
                annotations.append((repclass, subtype, overlap_pcrt))

    if len(annotations) == 0:
        return ('Simple Region', 'None', 0)

    sorted_annotations = sorted(annotations, key=lambda x: x[1], reverse=True)

    return sorted_annotations[0]


def main():
    remote_dir = '/Users/jiadonglin/Prostate/somatic_svs/'

    # separate_sniffles_somatic(remote_dir, 'T2T', 'winnowmap')
    add_haptag(f'{remote_dir}/Sniffles')

    annot_severus(f'{remote_dir}/severus')


    # prepare_pangenie_filter_files(remote_dir, 'GRCh38')

    # get_pop_filtered_somatics(remote_dir, 'GRCh38')
    # merge_somatics_each_sample(remote_dir, 'GRCh38')

    # annot_somatic_trans(remote_dir)
    # extract_somatic_ins(remote_dir)

    # ins_repeat_annot('/data/home/jdlin/Prostate/VCaP/ONT/sniffles/RM_INS_results',
    #                  'VCaP.somatic.ins.fasta.out.tab', 'VCaP.ins.RM.summary.txt')
    # vcap_pop_filtered_svs('/data/home/jdlin/Prostate/VCaP/ONT/sniffles', 'VCaP.minimap2.sniffles.PanGenie.merged.vcf')
    # extract_somatic_ins_seq('/data/home/jdlin/Prostate/VCaP/ONT/sniffles', '/data/home/jdlin/Prostate/VCaP/ONT/sniffles/VCaP.PanGenie.filtered.vcf', 'VCaP')

if __name__ == '__main__':
    main()
