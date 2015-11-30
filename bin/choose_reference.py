#!/usr/bin/env python

# a tool to select the best reference genome given a set of reads

import argparse
import os
import random
import sys

import bias

# dependencies
# * seqtk
# * bwa

BAM_TO_SAM="samtools view -h %s"
BEDTOOLS="bedtools"
BOWTIE_PATH="bowtie2"
BWA_PATH="/mnt/work/reference-bias/tools/bwa-0.7.12/bwa"
#BWA_PATH="bwa"
MAUVE_PATH="progressiveMauve"
SAMTOOLS="samtools"

def run( cmd ):
  '''
    run a system command
  '''
  bias.log_stderr( cmd )
  os.system( cmd )

def index( mapper, fasta ):
  if mapper == 'bwa':
    run( '%s index %s' % ( BWA_PATH, fasta ) )
  elif mapper == 'bowtie2':
    run( '%s-build %s %s-bt2' % ( BOWTIE_PATH, fasta, fasta ) )

def align( mapper, fasta, fastq, sam_out ):
  if mapper == 'bwa':
    run( '%s mem -t 8 %s %s > %s' % ( BWA_PATH, fasta, fastq, sam_out ) )
  elif mapper == 'bowtie2':
    run( '%s --local -p 16 -x %s-bt2 -U %s --quiet -S %s' % ( BOWTIE_PATH, fasta, fastq, sam_out ) ) # -t adds time

def sample_fastq( src, target, count ):
  run( 'seqtk sample {0} {1} > {2}'.format( src, count, target ) )

def unmapped_count( fasta, fastq, mapper='bwa', skipindex=False ):
  sam_tmp = 'tmp{0}.sam'.format( random.randint(1, 1e6) )
  if not skipindex:
    index( mapper, fasta ) 
  align( mapper, fasta, fastq, sam_tmp )
  stats_tmp = 'tmp{0}.stat'.format( random.randint(1, 1e6) )
  run( 'samtools flagstat {0} > {1}'.format( sam_tmp, stats_tmp ) )
  try:
    for line in open( stats_tmp, 'r' ):
      fields = line.strip().split()
      if len(fields) > 3 and fields[3] == 'mapped':
        return int( fields[0] )
  finally:
    run( 'rm {0} {1}'.format( sam_tmp, stats_tmp ) )

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Choose reference')
  parser.add_argument('--fastq', help='set of reads')
  parser.add_argument('--sample', type=int, default=10000, help='number of reads to sample')
  parser.add_argument('--skipindex', action='store_true', default=False, help='skip fasta indexing step')
  parser.add_argument('references', nargs='+', help='reference fasta file(s)')
  args = parser.parse_args()

  # sample fastq
  fastq_tmp = 'tmp{0}.fq'.format( random.randint(1, 1e6) )
  sample_fastq( args.fastq, fastq_tmp, args.sample )

  # assess with unmapped 
  results = {}
  bias.log_stderr( '{0} reference sequences'.format( len(args.references) ) )
  for reference in args.references:
    results[reference] = unmapped_count( reference, fastq_tmp, skipindex=args.skipindex )

  for k in sorted( results, key=results.get ):
    sys.stdout.write( '{0}\t{1}\n'.format( k, results[k] ) )

  run( 'rm {0}'.format( fastq_tmp ) )
