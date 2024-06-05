#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/6/2

'''


def select_nodes_from_chroms(paf_file, out_csv, chroms):

    csv_writer = open(out_csv, 'w')

    qname_dict = {}
    for line in open(paf_file, 'r'):
        entries = line.strip().split('\t')
        qname, chrom, ref_start, ref_end = entries[0], entries[5], entries[7], entries[8]
        if chrom in chroms:
            if qname in qname_dict:
                qname_dict[qname].append((chrom, ref_start, ref_end))
            else:
                qname_dict[qname] = [(chrom, ref_start, ref_end)]

    for qname, bp_info in qname_dict.items():
        bps = ''
        for bp in bp_info:
            bps += f'{bp[0]}:{bp[1]}-{bp[2]};'
        print(f'{qname}\t{len(bp_info)}\t{bps[:-1]}', file=csv_writer)


def add_align_info_to_nodes(paf_file, out_csv):
    csv_writer = open(out_csv, 'w')
    qname_dict = {}
    for line in open(paf_file, 'r'):
        entries = line.strip().split('\t')
        qname, chrom, ref_start, ref_end = entries[0], entries[5], entries[7], entries[8]
        if qname in qname_dict:
            qname_dict[qname].append((chrom, ref_start, ref_end))
        else:
            qname_dict[qname] = [(chrom, ref_start, ref_end)]

    for qname, bp_info in qname_dict.items():
        bps = ''
        chr_set = set()
        for bp in bp_info:
            bps += f'{bp[0]}:{bp[1]}-{bp[2]};'
            chr_set.add(bp[0])
        chr_str = ';'.join(list(chr_set))
        print(f'{qname}\t{chr_str}\t{len(bp_info)}\t{bps[:-1]}', file=csv_writer)

def main():

    paf_file = '/data/home/jdlin/Prostate/E/assembly/wgs_hifiasm_utg/BY-B-phased_SNVIndel_by_ALL_hap1.bp.p_utg.paf'
    csv_out = '/data/home/jdlin/Prostate/E/assembly/wgs_hifiasm_utg/utigs_at_chr12_chr17.csv'
    chroms = ['chr12', 'chr17']
    # select_nodes_from_chroms(paf_file, csv_out, chroms)

    add_align_info_to_nodes(paf_file, '/data/home/jdlin/Prostate/E/assembly/wgs_hifiasm_utg/utigs_align_info.csv')

if __name__ == '__main__':
    main()