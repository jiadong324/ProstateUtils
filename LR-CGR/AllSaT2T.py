#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/7/3

'''

import gzip
import json
import sys,os
import subprocess
import numpy as np

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from Helpers.Functions import *


def collect_signatures(outdir, prefix, bam_file):

    rmsk = pysam.Tabixfile('/data/home/jdlin/reference/T2T/hs1.repeatMasker.sorted.bed.gz', 'r')
    censat = pysam.Tabixfile('/data/home/jdlin/reference/T2T/chm13v2.0_censat_v2.0.bed.gz', 'r')


    aln_file = pysam.AlignmentFile(bam_file)
    aligns = aln_file.fetch()
    read_nr = 0

    sa_dict = {}
    while True:
        try:
            current_alignment = next(aligns)

            if current_alignment.is_unmapped or current_alignment.is_secondary:
                continue

            read_nr += 1
            # if read_nr % 10000 == 0:
            #     print("\tProcessed read {0}".format(read_nr))

            supplementary_alignments = retrieve_other_alignments(current_alignment, aln_file)
            good_suppl_alns = [aln for aln in supplementary_alignments if not aln.is_unmapped]
            sorted_sigs = analyze_trans_segments(aln_file, current_alignment, good_suppl_alns)

            if len(sorted_sigs) >= 2:
                sig_out, sig_count = [], 0

                for sig in sorted_sigs:
                    sig_count += 1
                    sig_chrom, sig_start, sig_end = aln_file.get_reference_name(sig['ref_id']), sig['ref_start'], sig['ref_end']
                    repclass, repsubtype, pcrt = rmsk_annotation(sig_chrom, sig_start, sig_end, 0, rmsk)
                    cen, pcrt = censat_annotation(sig_chrom, sig_start, sig_end, 0, censat)

                    ori = '+' if not sig['is_reverse'] else '-'
                    sig_out.append((sig_chrom, sig_start, sig_end, repclass, cen, ori, current_alignment.mapping_quality, sig['q_start'], sig['q_end']))

                sa_dict[current_alignment.query_name] = f'{current_alignment.query_name}\t{len(sorted_sigs)}\t{sig_out[:-1]}'
                sa_dict[current_alignment.query_name] = sig_out

        except StopIteration:
            break

    sa_bp_out = open(f'{outdir}/{prefix}.bed', 'w')
    for read_name, sa_info in sa_dict.items():
        for sa in sa_info:
            print(f'{sa[0]}\t{sa[1]}\t{sa[2]}\t{sa[3]}\t{sa[4]}\t{sa[5]}\t{sa[6]}\t{read_name}', file=sa_bp_out)
    sa_bp_out.close()

    hout = open(os.path.join(outdir, f"{prefix}.sorted.bed"), 'w')
    subprocess.check_call(["sort", "-k1,1", "-k2,2n", "-k3,3n", f'{outdir}/{prefix}.bed'], stdout=hout)
    hout.close()

    gout = open(os.path.join(outdir, f"{prefix}.sorted.bed.gz"), 'w')
    subprocess.check_call(["bgzip", "-f", "-c", os.path.join(outdir, f"{prefix}.sorted.bed")], stdout=gout)
    hout.close()

    subprocess.check_call(["tabix", "-p", "bed", os.path.join(outdir, f"{prefix}.sorted.bed.gz")])
    subprocess.check_call(["rm", "-rf", os.path.join(outdir, f"{prefix}.bed")])
    subprocess.check_call(["rm", "-rf", os.path.join(outdir, f"{prefix}.sorted.bed")])

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

        # if i > 0:
        #     prev_aln = alignments[i - 1]
        #     if prev_aln.reference_id == alignment.reference_id:
        #         if abs(alignment.reference_start - prev_aln.reference_start) > 100000:
        #             ref_dist = abs(alignment.reference_start - prev_aln.reference_start)

        new_alignment_dict = {'q_start': q_start,
                              'q_end': q_end,
                              'ref_id': alignment.reference_id,
                              'ref_start': alignment.reference_start,
                              'ref_end': alignment.reference_end,
                              'is_reverse': alignment.is_reverse}

        alignment_list.append(new_alignment_dict)

    sorted_alignment_list = sorted(alignment_list, key=lambda aln: (aln['q_start'], aln['q_end']))

    # if ref_dist < 100000 and len(chrom_set) == 1:
    #     return []

    return sorted_alignment_list

def detect_somatic_connections(workdir, bam_path, max_bpshift, min_sr, min_mapq):
    control_junction = pysam.Tabixfile(f'{workdir}/control_trans_junctions.sorted.bed.gz')
    tumor_junction_file = f'{workdir}/tumor_trans_junctions.sorted.bed.gz'
    aln_file = pysam.AlignmentFile(bam_path)

    somatic_junc_reads = {}
    somatic_junc_bps = {chrom: [] for chrom in VALID_CHROMS}
    for line in gzip.open(tumor_junction_file, 'rt'):
        chrom, ref_start, ref_end, rmsk, censat, ori, mapq, qname = line.strip().split('\t')
        # ref_pos = int(ref_end) if ori == '+' else int(ref_start)

        if chrom not in VALID_CHROMS:
            continue

        if censat == 'None_censat' and int(mapq) < min_mapq:
            continue

        if is_control_junction(chrom, int(ref_start), int(ref_end), control_junction, min_mapq):
            continue

        if qname in somatic_junc_reads:
            somatic_junc_reads[qname].append({'chrom': chrom, 'rstart': ref_start, 'rend': ref_end, 'rmsk': rmsk, 'censat': censat, 'ori': ori, 'mapq': mapq})
        else:
            somatic_junc_reads[qname] = [{'chrom': chrom, 'rstart': ref_start, 'rend': ref_end, 'rmsk': rmsk, 'censat': censat, 'ori': ori, 'mapq': mapq}]


        somatic_junc_bps[chrom].append({'chrom':chrom, 'rstart': int(ref_start), 'rend': int(ref_end), 'ori': ori,
                                        'mapq': mapq, 'qname': qname, 'rmsk': rmsk, 'censat': censat})

    filtered_somatic_junc = open(f'{workdir}/somatic_trans_junctions_sreads.bed', 'w')

    for qname, bp_list in somatic_junc_reads.items():
        bp_out = []

        for bp in bp_list:
            bp_out.append('{0}:{1}-{2},{3},{4},{5},{6}'.format(bp['chrom'], bp['rstart'], bp['rend'],  bp['ori'], bp['mapq'], bp['rmsk'], bp['censat']))
        bp_out_str = '\t'.join(bp_out)

        print(f'{qname}\t{len(bp_list)}\t{bp_out_str}', file=filtered_somatic_junc)

    with open(f'{workdir}/somatic_trans_junctions_sreads.json', 'w') as f:
        json.dump(somatic_junc_reads, f, indent=2)


    all_clusters = []
    for chrom, bp_list in somatic_junc_bps.items():
        this_clusters = group_trans_bp(bp_list, max_bpshift, min_sr)
        all_clusters.extend(this_clusters)

    print(f'\t# Detected split-read breakpoints {len(all_clusters)}')
    cluster_writer = open(f'{workdir}/somatic_trans_junction_clusters.tsv', 'w')

    link_info_dict = connect_clusters(all_clusters, aln_file, cluster_writer)

    link_writer = open(f'{workdir}/somatic_trans_junctions_links.tsv', 'w')

    link_counter = 0
    for link_id, link_info_list in link_info_dict.items():
        if len(link_info_list) > min_sr:
            link_counter += 1
            link_pos, link_rmsk, link_censa = ';'.join(link_info_list[0]['pos_info']), ';'.join(link_info_list[0]['rmsk']), ';'.join(link_info_list[0]['censat'])
            link_cov = ';'.join([str(ele) for ele in link_info_list[0]['pos_cov']])
            link_mapq = ';'.join([str(ele) for ele in link_info_list[0]['pos_mapq']])

            link_qnames = ';'.join([ele['qname'] for ele in link_info_list])
            link_sreads = ';'.join([ele['ref'] for ele in link_info_list])

            brpk_pos_af = []
            for brpk_cov in link_info_list[0]['pos_cov']:
                af = round(len(link_info_list) / brpk_cov, 2)
                brpk_pos_af.append(str(af))

            brpk_pos_af_str = ';'.join(brpk_pos_af)
            print(f'{link_id}\t{len(link_info_list)}\t{link_pos}\t{link_mapq}\t{link_cov}\t{brpk_pos_af_str}\t{link_rmsk}\t{link_censa}\t{link_qnames}\t{link_sreads}', file=link_writer)

    print(f'\t# Breakpoint links: {link_counter}')

def is_control_junction(chrom, t_start, t_end , control_junctions, min_mapq):

    for record_line in control_junctions.fetch(chrom, t_start, t_end):
        entries = record_line.strip().split('\t')
        c_start, c_end, mapq = int(entries[1]), int(entries[2]), int(entries[6])
        # if reciprocal_overlap(c_start, c_end, t_start, t_end) >= 0.8:
        if mapq < min_mapq:
            continue
        if (c_start - 100 < t_start and t_start < c_start + 100) or (c_end - 100 < t_end and t_end < c_end + 100):
            return True

    return False

def connect_clusters(all_bp_clusters, aln_file, cluster_writer):

    qname_clusters = {}
    for idx, cluster in enumerate(all_bp_clusters):
        cluster_chrom, cluster_brpk, mapq_median, bp_clusters = cluster
        cluster_qname_info = []

        try:
            cluster_brpk_cov = aln_file.count(cluster_chrom, cluster_brpk, cluster_brpk + 1)
            if cluster_brpk_cov == 0:
                cluster_brpk_cov = aln_file.count(cluster_chrom, cluster_brpk - 1, cluster_brpk)

            for bp in bp_clusters:
                qname = bp['qname']
                bp_info = '{0}:{1}-{2},{3},{4}'.format(bp['chrom'], bp['rstart'], bp['rend'], bp['ori'], bp['qname'])
                cluster_qname_info.append(bp_info)
                if qname in qname_clusters:
                    qname_clusters[qname].append((idx, cluster_chrom, cluster_brpk, mapq_median, cluster_brpk_cov, bp_clusters))
                else:
                    qname_clusters[qname] = [(idx, cluster_chrom, cluster_brpk, mapq_median, cluster_brpk_cov, bp_clusters)]
            cluster_qname_info_out = ';'.join(cluster_qname_info)
            print(f'{idx}\t{len(bp_clusters)}\t{cluster_chrom}\t{cluster_brpk}\t{mapq_median}\t{cluster_brpk_cov}\t{cluster_qname_info_out}', file=cluster_writer)
        except Exception:
            print(f'Coverage calculate error at {cluster_chrom}:{cluster_brpk}')

    link_info = {}
    for qname, clusters in qname_clusters.items():
        if contain_valid_connection(clusters):
            link_id = ''
            pos_sread = ''
            pos_info = []
            pos_cov = []
            pos_mapq = []
            pos_rmsk = []
            pos_censa = []

            for (cluster_id, cluster_chr, cluster_brpk, mapq_median, cluster_brpk_cov, cluster_sreads) in clusters:
                link_id += f'{cluster_id}-'
                # this_cluster_start, this_cluster_end = cluster_bps[0]['rstart'], cluster_bps[0]['rend']
                this_cluster_pos = f'{cluster_chr}:{cluster_brpk}'
                pos_cov.append(cluster_brpk_cov)
                pos_mapq.append(mapq_median)
                cluster_rmsk = set()
                cluster_censa = set()
                for ssv_read in cluster_sreads:
                    if ssv_read['qname'] == qname:
                        ori = 'F' if ssv_read['ori'] == '+' else 'R'
                        pos_sread = '{0}-{1}-{2}-{3}-{4}'.format(ssv_read['chrom'], ssv_read['rstart'], ssv_read['rend'], ori, ssv_read['mapq'])
                    cluster_rmsk.add(ssv_read['rmsk'])
                    cluster_censa.add(ssv_read['censat'])

                pos_info.append(this_cluster_pos)
                pos_rmsk.append(','.join(list(cluster_rmsk)))
                pos_censa.append(','.join(list(cluster_censa)))

            if link_id[:-1] in link_info:
                link_info[link_id[:-1]].append({'pos_info': pos_info, 'pos_mapq': pos_mapq,  'pos_cov': pos_cov,'rmsk': pos_rmsk, 'censat': pos_censa, 'qname': qname, 'ref': pos_sread})
            else:
                link_info[link_id[:-1]] = [{'pos_info': pos_info, 'pos_mapq': pos_mapq, 'pos_cov': pos_cov, 'rmsk': pos_rmsk, 'censat': pos_censa, 'qname': qname, 'ref': pos_sread}]

    return link_info

def group_trans_bp(bp_list, max_bpshift, min_sr):
    rstart_sorted = sorted(bp_list, key=lambda x: x['rstart'])
    clusters = []
    this_cluster = [rstart_sorted[0]]

    cluster_pos_tracker = {}
    for i in range(1, len(rstart_sorted)):
        prev_bp_start, curr_bp_start = this_cluster[-1]['rstart'], rstart_sorted[i]['rstart']
        if prev_bp_start - max_bpshift > curr_bp_start or prev_bp_start + max_bpshift < curr_bp_start:
            if len(this_cluster) > min_sr:

                this_cluster_chrom, this_cluster_start, mapq_median = get_left_brpk_from_cluster(this_cluster)
                clusters.append([this_cluster_chrom, this_cluster_start, mapq_median, this_cluster])

                cluster_pos_tracker[f'{this_cluster_chrom}:{this_cluster_start}'] = 1

            this_cluster = []

        this_cluster.append(rstart_sorted[i])

    rend_sorted = sorted(bp_list, key=lambda x: x['rend'])

    this_cluster = [rend_sorted[0]]
    for i in range(1, len(rend_sorted)):
        prev_bp_start, prev_bp_end, curr_bp_start, curr_bp_end = this_cluster[-1]['rstart'], this_cluster[-1]['rend'], rend_sorted[i]['rstart'], rend_sorted[i]['rend']
        if prev_bp_end - max_bpshift > curr_bp_end or prev_bp_end + max_bpshift < curr_bp_end:

            if len(this_cluster) > min_sr and not is_cluster_brpk_duplicated(this_cluster, cluster_pos_tracker):
                cluster_chrom, cluster_brpk, cluster_mapq = get_right_brpk_from_cluster(this_cluster)
                clusters.append([cluster_chrom, cluster_brpk, cluster_mapq, this_cluster])

            this_cluster = []
        this_cluster.append(rend_sorted[i])

    return clusters

def is_cluster_brpk_duplicated(this_cluster, cluster_pos_tracker):
    cluster_chrom, cluster_brpk, cluster_mapq = get_right_brpk_from_cluster(this_cluster)
    return True if f'{cluster_chrom}:{cluster_brpk}' in cluster_pos_tracker else False

def contain_valid_connection(clusters):
    if len(clusters) == 1:
        return False

    sorted_clusters = sorted(clusters, key=lambda x:x[2])
    cluster_span = sorted_clusters[-1][2] - sorted_clusters[0][2]
    cluster_chroms = set()
    for cluster in clusters:
        cluster_chroms.add(cluster[1])

    if len(cluster_chroms) > 1:
        return True

    if len(cluster_chroms) == 1 and cluster_span > 100000:
        return True

    return False

def get_left_brpk_from_cluster(cluster):
    bp_info_dict = {}
    mapqs = []
    for bp in cluster:
        pos1 = bp['rstart']
        mapqs.append(int(bp['mapq']))
        if pos1 in bp_info_dict:
            bp_info_dict[pos1].append(bp)
        else:
            bp_info_dict[pos1] = [bp]
    sorted_bp_info_dict = sorted(bp_info_dict.items(), key=lambda x: len(x[1]), reverse=True)

    return cluster[0]['chrom'], sorted_bp_info_dict[0][0], int(np.median(mapqs))

def get_right_brpk_from_cluster(cluster):
    bp_info_dict = {}
    mapqs = []
    for bp in cluster:
        pos2 = bp['rend']
        mapqs.append(int(bp['mapq']))

        if pos2 in bp_info_dict:
            bp_info_dict[pos2].append(bp)
        else:
            bp_info_dict[pos2] = [bp]
    sorted_bp_info_dict = sorted(bp_info_dict.items(), key=lambda x: len(x[1]), reverse=True)

    return cluster[0]['chrom'], sorted_bp_info_dict[0][0], int(np.median(mapqs))


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
    # sample = 'B'
    aligner = 'minimap2'
    control_junction_minq = 20
    somatic_bp_minsupp = 3
    somatic_bp_maxshift = 10

    # for sample in SAMPLES:
    #     workdir = f'/data/home/jdlin/Prostate/{sample}/somatic_trans/T2T/{aligner}/all_sa'
    #
    #     if not os.path.exists(workdir):
    #         os.mkdir(workdir)
    #
    #     bamdir = f'/data/DATA/PRAD-CN/Li_PRAD/{sample}'
    #
    #     for tissue in ['T', 'B']:
    #         prefix = 'control_trans_junctions' if tissue == 'B' else 'tumor_trans_junctions'
    #         bam_suffix = 'minimap2.haptag.sorted.bam' if tissue == 'T' else 'minimap2.sorted.bam'
    #
    #         if aligner == 'winnowmap':
    #             bam_suffix = 'winnowmap.sorted.bam'
    #         print(f'Collect junctions {sample} {prefix}')
    #         collect_signatures(workdir, prefix, f'{bamdir}/{sample}-{tissue}.T2T.{bam_suffix}')
    #
    #     print(f'Detect somatic junctions from: {sample} ...')
    #     bam_suffix = 'winnowmap.sorted.bam' if aligner == 'winnowmap' else 'minimap2.haptag.sorted.bam'
    #     detect_somatic_connections(workdir, f'{bamdir}/{sample}-T.T2T.{bam_suffix}', somatic_bp_maxshift, somatic_bp_minsupp, control_junction_minq)

    bam_file = sys.argv[1]
    prefix = sys.argv[2]
    outdir = sys.argv[3]

    # collect_signatures(outdir, prefix, bam_file)
    detect_somatic_connections(outdir, bam_file, somatic_bp_maxshift,
                               somatic_bp_minsupp, control_junction_minq)

if __name__ == '__main__':
    main()