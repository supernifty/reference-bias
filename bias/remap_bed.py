#!/usr/bin/env python

# given an xmfa, remap a bam file coordinates
# usage remap_bed --xmfa xmfa --origin 0 --target 1 srcbamfile > targetsamfile
# writes coverage stats to stderr

import argparse
import sys

import bias

def remap_bed(xmfa, origin, target, remap_cigar, bed, chromosome, out_fh):
  # measure coverage
  sys.stderr.write("building coverage map for bed remap...\n")
  #print args.origin, dest
  mauve_map = bias.MauveMap(open(xmfa, 'r'), src_strand=origin, target_strand=target)
  sys.stderr.write( "{} positions mapped\n".format(len(mauve_map.coverage)))

  sys.stderr.write("analyzing bed {0}...\n".format( bed ))
  mapped = 0
  total = 0
  for line in open(bed, 'r'):
      fields = line.strip('\n').split('\t')
      if fields[0] == chromosome:
          gene = fields[3]
          for pos in range(int(fields[1]), int(fields[2])):
              if pos in mauve_map.coverage: # we can map it
                  mapped += 1
                  new_pos = mauve_map.coverage[pos]
                  out_fh.write('{}\t{}\t{}\t{}\n'.format(fields[0], new_pos, new_pos+1, gene))
              total += 1

  sys.stderr.write("analyzing bed {0}: mapped {1} positions out of {2}\n".format(bed, mapped, total))

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Remap a bam file using XMFA')
  parser.add_argument('bed', help='bed file to analyze')
  parser.add_argument('--xmfa', required=True, help='xmfa file')
  parser.add_argument('--origin', type=int, default=1, help='index of fasta that bam was mapped to')
  parser.add_argument('--target', type=int, default=2, help='index of fasta that bam should be mapped to')
  parser.add_argument('--remap_cigar', required=False, default=False, help='attempt to remap the cigar string')
  parser.add_argument('--chromosome', help='filter bed on chromosome')

  args = parser.parse_args()
  remap_bed( args.xmfa, args.origin, args.target, args.remap_cigar, args.bed, args.chromosome, sys.stdout)
 
