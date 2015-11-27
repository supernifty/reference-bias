#!/usr/bin/env python

# given input from compare_bams, generate a sam file with matching mismatches based on specified criteria
# e.g. python extract_mismatched_reads.py --min_distance 50 --mappable True ../../data/e-coli-mg1655_snp_ins_del_50_x1_bwa.bam < ../../data/compare_bams3.log > mappable_test.sam

import argparse
import sys

import bias
#import config

def extract_mismatched_reads( min_distance, max_distance, mappable, bam, in_fh, out_fh, bam_to_sam ):
  bias.log_stderr( 'Filtering based on min distance %i, max_distance %i, mappable %s' % ( min_distance, max_distance, mappable ) )

  allowed_tags = set()
  candidates = 0
  for line in in_fh:
    fields = line.strip().split(',')
    if len(fields) == 5 and fields[0].isdigit():
      candidates += 1
      distance = abs(int(fields[3]))
      line_mappable = int(fields[4]) == 0
      if distance >= min_distance and distance <= max_distance:
        if line_mappable and mappable or not line_mappable and not mappable:
          allowed_tags.add( fields[2] )
  bias.log_stderr( '%i allowed tags from %i possibles' % ( len(allowed_tags), candidates ) )

  bias.SamFilter( sam_fh=bias.BamReaderExternal( bam_to_sam, bam ), target_fh=out_fh, allowed_tags=allowed_tags, log=bias.log_stderr )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Find mismatched reads matching criteria')
  parser.add_argument('bam', help='bam files to extract reads')
  parser.add_argument('--min_distance', dest='min_distance', type=int, default=0, help='only consider alignments with at least this distance from the true location')
  parser.add_argument('--max_distance', dest='max_distance', type=int, default=1e12, help='only consider alignments with at most this distance from the true location')
  parser.add_argument('--mappable', dest='mappable', type=bool, default=False, help='only consider alignments that have a mappable location')

  args = parser.parse_args()

  extract_mismatched_reads( args.min_distance, args.max_distance, args.mappable, args.bam, sys.stdin, sys.stdout )
