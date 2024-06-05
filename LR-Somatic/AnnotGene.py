#!/usr/bin/env python

# encoding: utf-8

'''
@author: Jiadong Lin

@contact: jiadong324@gmail.com

@time: 2022/12/7

'''

import os
import pysam

from Helpers.Constants import *

def bp_at_cosmic_genes(workdir):
    cosmic_genes = {}
    with open(REMOTECOSMIC, 'r') as f:
        next(f)
        for line in f:
            entries = line.strip().split('\t')
            gene, coord = entries[0], entries[3]

            try:
                chr, start, end = 'chr' + coord.split(':')[0], int(coord.split(':')[1].split('-')[0]), int(coord.split(':')[1].split('-')[1])
                if chr in cosmic_genes:
                    cosmic_genes[chr][gene] = (start, end)
                else:
                    cosmic_genes[chr] = {}
                    cosmic_genes[chr][gene] = (start, end)
            except Exception:
                print(line)

    # for sample in SAMPLES:



def main():
    workdir = '/data/home/jdlin/Prostate'

    bp_at_cosmic_genes(workdir)


if __name__ == '__main__':

    main()