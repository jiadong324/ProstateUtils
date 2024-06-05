#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/5/31

'''
import sys,os

import pysam
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Functions import *


def rename_normal_haps(hap1_fa, hap2_fa, outdir):
    h1_tagged = []
    h2_tagged = []
    merged_haps = []



    for h1_rec in list(SeqIO.parse(hap1_fa, "fasta")):
        rec_name, rec_seq = h1_rec.id, h1_rec.seq
        new_rec = SeqRecord(rec_seq, id=f'H1-{rec_name}', description=f'H1-{rec_name}', name=f'H1-{rec_name}')
        merged_haps.append(new_rec)
        h1_tagged.append(new_rec)

    for h2_rec in list(SeqIO.parse(hap2_fa, "fasta")):
        rec_name, rec_seq = h2_rec.id, h2_rec.seq
        new_rec = SeqRecord(rec_seq, id=f'H2-{rec_name}', description=f'H2-{rec_name}', name=f'H2-{rec_name}')
        merged_haps.append(new_rec)
        h2_tagged.append(new_rec)

    # merged_fa = f'{outdir}/normal.h1-h2.ragtag.scaffold.fasta'

    # with open(out_fa, 'w') as out:
    #     SeqIO.write(merged_haps, out, 'fasta')


    with open(f'{outdir}/h1/ragtag.scaffold.tagged.h1.fasta', 'w') as fout:
        SeqIO.write(h1_tagged, fout, 'fasta')


    with open(f'{outdir}/h2/ragtag.scaffold.tagged.h2.fasta', 'w') as fout:
        SeqIO.write(h2_tagged, fout, 'fasta')

def separate_ref_aln_bam(input_bam, outdir, out_prefix):
    aln_file = pysam.AlignmentFile(input_bam)
    aligns = aln_file.fetch()

    h1_bam = pysam.AlignmentFile(f'{outdir}/{out_prefix}.H1.bam', 'wb', template=aln_file)
    h2_bam = pysam.AlignmentFile(f'{outdir}/{out_prefix}.H2.bam', 'wb', template=aln_file)

    while True:
        try:
            current_alignment = next(aligns)
            ref_id = current_alignment.reference_id
            ref_name = aln_file.get_reference_name(ref_id)

            if ref_name.split('-')[0] == 'H1':
                h1_bam.write(current_alignment)

            if ref_name.split('-')[0] == 'H2':
                h2_bam.write(current_alignment)

        except StopIteration:
            break

def select_abnormal_reads(workdir, bam_path, read_file, chrom, start, end):

    reads_list = []
    for line in open(read_file, 'r'):
        reads_list.extend(line.strip().split(','))

    aln_file = pysam.AlignmentFile(bam_path)

    aligns = aln_file.fetch(chrom, start, end)
    abnomal_records = []

    while True:
        try:
            current_alignment = next(aligns)

            if current_alignment.query_name in reads_list:
                abnomal_records.append(SeqRecord(Seq(current_alignment.query_sequence), id=current_alignment.query_name, description=current_alignment.query_name))

        except StopIteration:
            break

    with open(f'{workdir}/brpks_reads.fasta', 'w') as output:
        SeqIO.write(abnomal_records, output, "fasta")


# def reads_span_mutations(snp_vcf, mpileup_file):


def analyze_trans_segments(bam_file, primary, supplementaries):
    alignments = [primary] + supplementaries
    alignment_list = []
    chrom_set = set()
    ref_dist = 0
    for i, alignment in enumerate(alignments):
        # correct query coordinates for reversely mapped reads
        if alignment.is_reverse:
            q_start = alignment.infer_read_length() - alignment.query_alignment_end
            q_end = alignment.infer_read_length() - alignment.query_alignment_start
        else:
            q_start = alignment.query_alignment_start
            q_end = alignment.query_alignment_end

        chrom = bam_file.get_reference_name(alignment.reference_id)
        if chrom not in VALID_CHROMS:
            continue

        chrom_set.add(chrom)

        if i > 0:
            prev_aln = alignments[i - 1]
            if prev_aln.reference_id == alignment.reference_id:
                if abs(alignment.reference_start - prev_aln.reference_start) > 100000:
                    ref_dist = abs(alignment.reference_start - prev_aln.reference_start)

        new_alignment_dict = {'q_start': q_start,
                              'q_end': q_end,
                              'ref_id': alignment.reference_id,
                              'ref_start': alignment.reference_start,
                              'ref_end': alignment.reference_end,
                              'is_reverse': alignment.is_reverse}

        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))

    if ref_dist < 100000 and len(chrom_set) == 1:
        return []

    return sorted_alignment_list

def main():
    bamdir = '/data/DATA/PRAD-CN/Li_PRAD'
    sample, tissue = 'E', 'T'
    workdir = f'/data/home/jdlin/Prostate/{sample}/assembly'

    for sample in SAMPLES:
        rename_normal_haps(f'/data/home/jdlin/Prostate/{sample}/assembly/map2blood_haps/ragtag/h1/ragtag.scaffold.fasta',
                            f'/data/home/jdlin/Prostate/{sample}/assembly/map2blood_haps/ragtag/h2/ragtag.scaffold.fasta',
                            f'/data/home/jdlin/Prostate/{sample}/assembly/map2blood_haps/ragtag')


    # separate_ref_aln_bam('/data/home/jdlin/Prostate/E/assembly/map2blood_haps/ragtag/E-T.h1h2.minimap2.sorted.bam',
    #                      '/data/home/jdlin/Prostate/E/assembly/map2blood_haps/ragtag', 'E.tumor-reads')

    # select_abnormal_reads(f'{workdir}/brpk_reads_29620809', f'{bamdir}/{sample}/{sample}-{tissue}.GRCh38.ngmlr.sorted.addrg.haptag.bam',
    #                       f'{workdir}/brpk_reads_29620809/abnormal_reads_at_29635809.txt', 'chr12', 29634809, 29636809)

if __name__ == '__main__':
    main()