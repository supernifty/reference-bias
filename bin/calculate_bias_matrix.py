#!/usr/bin/env python

# provide a list of read sets and reference genomes
# generate a bias matrix

# python calculate_bias_matrix target_file < config
# config is space separated lines of the form
# sra,fasta

import itertools
import sys

import bias

BAM_TO_SAM="samtools view -h %s"
BEDTOOLS="bedtools"
BOWTIE_PATH="bowtie2"
BWA_PATH="/mnt/work/reference-bias/tools/bwa-0.7.12/bwa"
#BWA_PATH="bwa"
MAUVE_PATH="progressiveMauve"
SAMTOOLS="samtools"

def bias_matrix( refs, log, out ):
  out.write( 'Donor\tReference\tSRA\tLow\tMid\tHigh\n' )
  for donor, reference in itertools.product( refs, repeat=2 ):
    log.write( 'Calculating donor {0} reference {1}\n'.format( donor, reference ) )
    calculator = bias.Calculator( BWA_PATH, BOWTIE_PATH, MAUVE_PATH, BAM_TO_SAM, log, log )
    low, mid, high = calculator.calculate( donor=donor[1], reference=reference[1], job=None, stage=None, tmpdir='./tmp', align='bwa', donorbam=None, donorsam=None, fastq=donor[0] )
    #low, mid, high = 1, 2, 3
    out.write( '{0}\t{1}\t{2}\t{3:.1f}\t{4:.1f}\t{5:.1f}\n'.format( donor[1], reference[1], donor[0], low, mid, high ) )

if __name__ == '__main__':
  refs = []
  for line in sys.stdin:
    if line.startswith('#'):
      continue
    sra, ref = [ x.strip() for x in line.strip().split(',') ]
    refs.append( (sra, ref) )
  target_file = sys.argv[1]
  m = bias_matrix( refs, sys.stderr, open( target_file, 'w' ) ) 
