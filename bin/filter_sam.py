#!/usr/bin/env python

import argparse

def filter_sam( out_fn, in_fn, chromosome):
  with open(out_fn, 'w') as donor_out:
    for line in open(in_fn, 'r'):
      if line.startswith("@SQ"):
        if "SN:{}\t".format(chromosome) in line:
          donor_out.write(line)
      elif line.startswith("@"):
        donor_out.write(line)
      else:
        fields = line.strip('\n').split('\t')
        if fields[2] == chromosome:
          donor_out.write(line)

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Filter specified chromosome')
  parser.add_argument('--target', help='target filename')
  parser.add_argument('--source', help='source filename')
  parser.add_argument('--chromosome', help='name of chromosome')
  args = parser.parse_args()
  filter_sam(args.target, args.source, args.chromosome)
