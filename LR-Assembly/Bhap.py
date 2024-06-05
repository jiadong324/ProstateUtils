#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2023/6/19

'''
import networkx as nx
import json


class Edge:
    def __init__(self, node1, node2):
        self.node1 = node1
        self.node2 = node2

    def __eq__(self, other):
        return (self.node1 == other.node1 and self.node2 == other.node2) or (self.node1 == other.node2 and self.node2 == other.node1)
        # return self.node1 == other.node2 and self.node2 == other.node1

    def __hash__(self):
        return hash(self.node1) ^ hash(self.node2)

    def tostring(self):
        return f'{self.node1}-{self.node2}'

    def get_node1(self):
        return self.node1

    def get_node2(self):
        return self.node2


def parse_blood_gfa(workdir, gfa_input):
    node_info_dict = {}
    for line in open(gfa_input, 'r'):
        entries = line.strip().split('\t')
        if entries[0] == 'S':
            node_id, node_chrom, node_len, node_seq = entries[1], entries[4].split(':')[2], int(entries[3].split(':')[2]), entries[2]
            node_info_dict[node_id] = (node_chrom, node_len, node_seq)

    with open(f'{workdir}/ragtag.scaffold.h1h2.nodes.json', 'w') as f:
        json.dump(node_info_dict, f)

    print('Loading normal haplotype graph ...')

    G = nx.Graph()
    edge_counter = 0
    for line in open(gfa_input, 'r'):
        entries = line.strip().split('\t')
        if entries[0] == 'L':
            edge_counter += 1
            node1, node1_ori, node2, node2_ori = entries[1], entries[2], entries[3], entries[4]
            node1_len, node2_len = entries[7].split(':')[2], entries[8].split(':')[2]
            node1_sn, node2_sn = node_info_dict[node1][0], node_info_dict[node2][0]
            node1_info = f'{node1_ori}-{node1_len}-{node1_sn}'
            node2_info = f'{node2_ori}-{node2_len}-{node2_sn}'
            G.add_node(node1, info=node1_info)
            G.add_node(node2, info=node2_info)
            G.add_edge(node1, node2)

    print(f'# Nodes: {len(node_info_dict)}\n# Edges: {edge_counter}')

    nhap_cc = open(f'{workdir}/ragtag.scaffold.h1h2.cc.tsv', 'w')
    G_cc = sorted(nx.connected_components(G), key=len, reverse=True)
    for idx, cc in enumerate(G_cc):
        cc_size = len(cc)
        nodes_str = ','.join(list(cc))
        ragtag_sns = set()
        for node in list(cc):
            ragtag_sns.add(node_info_dict[node][0])
        sns_str = ','.join(list(ragtag_sns))
        print(f'cc-{idx}\t{sns_str}\t{cc_size}\t{nodes_str}', file=nhap_cc)
    print(f'# CC: {len(G_cc)}')


def parse_tumor_augmented_gfa(workdir, tumor_gaf, blood_gfa_nodes, blood_gfa, min_sup, max_sa, min_mapq):

    print('Loading normal haplotype graph ...')

    node_info_dict = json.load(open(blood_gfa_nodes))

    G = nx.Graph()
    edge_counter = 0
    for line in open(blood_gfa, 'r'):
        entries = line.strip().split('\t')
        if entries[0] == 'L':
            node1, node1_ori, node2, node2_ori = entries[1], entries[2], entries[3], entries[4]
            node1_len, node2_len = entries[7].split(':')[2], entries[8].split(':')[2]
            node1_sn, node2_sn = node_info_dict[node1][0], node_info_dict[node2][0]
            node1_info = f'{node1_ori}-{node1_len}-{node1_sn}'
            node2_info = f'{node2_ori}-{node2_len}-{node2_sn}'
            G.add_node(node1, info=node1_info)
            G.add_node(node2, info=node2_info)
            G.add_edge(node1, node2)
            # G.add_edges_from((node1, node2, {'type': 'control', 'weight': 1}))
            edge_counter += 1

    print('Adding tumor edges ...')

    tumor_reads_aligns = []
    tumor_reads_sa = {}
    tumor_reads_sa_info = {}
    new_edge_dict = {}

    for line in open(tumor_gaf, 'r'):
        entries = line.strip().split('\t')
        qname, qlen, qstart, qend, qori = entries[0].split(' ')[0], int(entries[1]), int(entries[2]), int(entries[3]), entries[4]
        node_path, rstart, rend, mres, mapq = entries[5], int(entries[7]), int(entries[8]), int(entries[9]), int(entries[11])

        if mapq < min_mapq:
            continue

        qalign_pcrt = (qend - qstart) / qlen * 100
        nodes_on_path = process_node_path(node_path)
        tumor_reads_aligns.append((qalign_pcrt, len(nodes_on_path)))


        exori_path = '-'.join(nodes_on_path)

        if qname in tumor_reads_sa:
            tumor_reads_sa[qname].extend(nodes_on_path)
            tumor_reads_sa_info[qname].append({'path': exori_path, 'qstart': qstart, 'qend': qend, 'rstart': rstart, 'rend': rend})
        else:
            tumor_reads_sa[qname] = nodes_on_path
            tumor_reads_sa_info[qname] = [{'path': exori_path, 'qstart': qstart, 'qend': qend, 'rstart': rstart, 'rend': rend}]


    ## Examining novel SSV connections
    for qname, node_list in tumor_reads_sa.items():
        # Too much alignments
        if len(tumor_reads_sa_info[qname]) > max_sa:
            continue
        if not contain_ssv_signature(tumor_reads_sa_info[qname]):
            continue

        for i in range(1, len(node_list)):
            edge_obj = Edge(node_list[i - 1], node_list[i])
            if not G.has_edge(node_list[i - 1], node_list[i]):
                G.add_edge(node_list[i - 1], node_list[i])
                new_edge_dict[edge_obj] = [qname]
            elif edge_obj in new_edge_dict:
                new_edge_dict[edge_obj].append(qname)

    added_edge_num = 0
    altered_node_set = set()
    added_edge_info = {}

    ## Write novel edge and alignments to file
    novel_edges = open(f'{workdir}/tumor_augmented_edges.tsv', 'w')
    for edge, qnames in new_edge_dict.items():
        if len(qnames) > min_sup:
            for qname in qnames:
                for aln in tumor_reads_sa_info[qname]:
                    if aln['path'] in added_edge_info:
                        added_edge_info[aln['path']].append({'qname': qname, 'qstart': aln['qstart'], 'qend': aln['qend'], 'rstart': aln['rstart'], 'rend': aln['rend']})
                    else:
                        added_edge_info[aln['path']] = [{'qname': qname, 'qstart': aln['qstart'], 'qend': aln['qend'], 'rstart': aln['rstart'], 'rend': aln['rend']}]

            qname_out = ';'.join(qnames)
            node1, node2 = edge.get_node1(), edge.get_node2()
            altered_node_set.add(node1)
            altered_node_set.add(node2)

            node1_ichr, node2_ichr = node_info_dict[node1][0], node_info_dict[node2][0]
            added_edge_num += 1
            print(f'{edge.tostring()}\t{node1_ichr};{node2_ichr}\t{len(qnames)}\t{qname_out}', file=novel_edges)

    path_aln_brpks = group_tumor_path_aln(added_edge_info, 50, min_sup)

    with open(f'{workdir}/tumor_augmented_aln_info.json', 'w') as f:
        json.dump(added_edge_info, f, indent=2)

    with open(f'{workdir}/tumor_augmented_path_brpks.json', 'w') as f:
        json.dump(path_aln_brpks, f, indent=2)

    print(f'# Tumor edges added: {added_edge_num}')
    print(f'# Tumor altered nodes: {len(altered_node_set)}')


def build_augmented_graph(chroms, blood_nodes_info, blood_gfa, blood_ccs, new_edge, path_brpks):

    # node_info_dict = json.load(open(blood_nodes_info))
    #
    # G = nx.FindPath()
    # for line in open(blood_gfa, 'r'):
    #     entries = line.strip().split('\t')
    #     if entries[0] == 'L':
    #         node1, node1_ori, node2, node2_ori = entries[1], entries[2], entries[3], entries[4]
    #         node1_len, node2_len = entries[7].split(':')[2], entries[8].split(':')[2]
    #         node1_sn, node2_sn = node_info_dict[node1][0], node_info_dict[node2][0]
    #         node1_info = f'{node1_ori}-{node1_len}-{node1_sn}'
    #         node2_info = f'{node2_ori}-{node2_len}-{node2_sn}'
    #         G.add_node(node1, info=node1_info)
    #         G.add_node(node2, info=node2_info)
    #         G.add_edge(node1, node2)
    #
    # blood_cc_dict = {}
    # for line in open(blood_ccs, 'r'):
    #     entries = line.strip().split('\t')
    #     cc_id, cc_chroms, nodes = entries[0], entries[1], entries[3]
    #     chrom_set = set()
    #     for cc_chrom in cc_chroms.split(','):
    #         chrom = cc_chrom.split('_')[0].split('-')[1]
    #         chrom_set.add(chrom)
    #     chrom_key = ','.join(list(chrom_set))
    #     blood_cc_dict[chrom_key] = nodes

    path_brpks_dict = json.load(open(path_brpks))

    for path, brpks_list in path_brpks_dict.items():
        nodes_on_path = process_node_path(path)



def group_tumor_path_aln(path_aln_dict, bp_shift, min_sr):

    group_aln_info = {}
    for path_name, aln_list in path_aln_dict.items():
        sorted_aln_list = sorted(aln_list, key=lambda x: x['rstart'])
        new_nodes = []
        # clusters = []
        this_cluster = [sorted_aln_list[0]]
        for i in range(1, len(sorted_aln_list)):
            last_aln, this_aln = this_cluster[-1], sorted_aln_list[i]
            abs_ref_dist = abs(this_aln['rstart'] - last_aln['rstart'])
            if abs_ref_dist > bp_shift:
                if len(this_cluster) > min_sr:
                    this_cluster_pos = this_cluster[0]['rstart']
                    this_cluster_rnames = ','.join([ele['qname'] for ele in this_cluster])
                    rstart = ','.join([str(ele['rstart']) for ele in this_cluster])
                    rend = ','.join([str(ele['rend']) for ele in this_cluster])
                    new_nodes.append({'rstart': rstart, 'rend': rend, 'pseudo_node_id': f'{path_name}_{len(new_nodes)}',
                                      'pseudo_node_rpos': this_cluster_pos, 'pseudo_node_supp': len(this_cluster), 'pseudo_node_rnames': this_cluster_rnames})

                    this_cluster = []
                else:
                    this_cluster = []

            this_cluster.append(this_aln)

        if len(this_cluster) > min_sr:
            this_cluster_pos = this_cluster[0]['rstart']
            this_cluster_rnames = ';'.join([ele['qname'] for ele in this_cluster])
            rstart = ','.join([str(ele['rstart']) for ele in this_cluster])
            rend = ','.join([str(ele['rend']) for ele in this_cluster])
            new_nodes.append({'rstart': rstart, 'rend': rend, 'pseudo_node_id': f'{path_name}_{len(new_nodes)}',
                              'pseudo_node_rpos': this_cluster_pos, 'pseudo_node_supp':len(this_cluster), 'pseudo_node_rnames': this_cluster_rnames})

        if len(new_nodes) > 0:
            group_aln_info[path_name] = new_nodes

    return group_aln_info


def contain_ssv_signature(read_aligns):

    sorted_read_aligns = sorted(read_aligns, key=lambda x: x['qstart'])

    has_sig = False
    for i in range(1, len(sorted_read_aligns)):
        prev_aln, curr_aln = sorted_read_aligns[i - 1], sorted_read_aligns[i]
        dist_on_read = curr_aln['qstart'] - prev_aln['qend']
        dist_on_ref = curr_aln['rstart'] - prev_aln['rend']

        if prev_aln['path'] == curr_aln['path']:
            # No overlap on ref
            if dist_on_ref >= -10:
                deviation = dist_on_read - dist_on_ref
                if deviation >= 50:
                    has_sig = True
            else:
                if dist_on_ref <= -50:
                    has_sig = True
        else:
            has_sig = True

    return has_sig


def process_node_path(node_path):
    nodes = []
    node_char = [node_path[0]]

    for i in range(1, len(node_path)):
        if node_path[i] == '>' or node_path[i] == '<':
            this_node = ''.join(node_char[1:])
            nodes.append(this_node)
            node_char = []

        if i == len(node_path) - 1:
            node_char.append(node_path[-1])
            this_node = ''.join(node_char[1:])
            nodes.append(this_node)
        node_char.append(node_path[i])

    return nodes

def reciprocal_overlap(begin_a, end_a, begin_b, end_b):

    overlap = min(end_a, end_b) - max(begin_a, begin_b)
    return round(min([overlap / (end_a - begin_a), overlap / (end_b - begin_b)]), 3)



def main():
    sample = 'E'
    workdir = f'/data/home/jdlin/Prostate/{sample}/assembly/map2blood_haps/ragtag'

    max_sa = 5
    min_sr = 5
    min_mapq = 30

    # parse_blood_gfa(workdir, f'{workdir}/ragtag.scaffold.h1h2.gfa')
    parse_tumor_augmented_gfa(workdir, f'{workdir}/tumor.h1h2.mg.graphaligner.gaf',
                              f'{workdir}/ragtag.scaffold.h1h2.nodes.json', f'{workdir}/ragtag.scaffold.h1h2.gfa', min_sr, max_sa, min_mapq)

    # build_augmented_graph(['chr12', 'chr17'], f'{workdir}/ragtag.scaffold.h1h2.nodes.json', f'{workdir}/ragtag.scaffold.h1h2.gfa',
    #                       f'{workdir}/ragtag.scaffold.h1h2.cc.tsv', f'{workdir}/tumor_augmented_edges.tsv', f'{workdir}/tumor_augmented_path_brpks.json')

if __name__ == '__main__':
    main()

