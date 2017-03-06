#!/usr/bin/env python

import argparse
import collections
import sys

def assess_bias(chrom, bias, bed):
    sys.stderr.write('reading bias file...')

    unaffected = set()
    total = set()

    first = True

    for idx, line in enumerate(bias):
        if idx % 1000000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        if first:
            first = False
            continue
        fields = line.strip('\n').split(',')
        # pos, correct, total
        if int(fields[1]) > 0:
            unaffected.add(int(fields[0]))
        total.add(int(fields[0]))

    sys.stderr.write('reading bed file...')
    gene_unaffected = collections.defaultdict(int)
    gene_total = collections.defaultdict(int)

    for idx, line in enumerate(bed):
        if idx % 100000 == 0:
            sys.stderr.write("read {} lines...\n".format(idx))
        fields = line.strip('\n').split('\t')

        # chr, start, end, gene
        if fields[0] == chrom:
            for pos in range(int(fields[1]), int(fields[2])):
                if pos in unaffected:
                    gene_unaffected[fields[3]] += 1
                if pos in total:
                    gene_total[fields[3]] += 1

    # write result
    sys.stdout.write('{},{},{}\n'.format('bias', 'gene', 'count'))
    for gene in gene_total.keys():
        sys.stdout.write('{:.1f},{}\n'.format(100 - 100 * gene_unaffected[gene] / gene_total[gene], gene, gene_total[gene]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Assess bias by feature')
    parser.add_argument('--bias', help='generated bias file')
    parser.add_argument('--bed', help='features')
    parser.add_argument('--chromosome', help='chromosome')
    args = parser.parse_args()
    assess_bias(args.chromosome, open(args.bias, 'r'), open(args.bed, 'r'))
