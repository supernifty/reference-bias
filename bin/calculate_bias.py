#!/usr/bin/env python
#####################################################
# the entire process of aligning and calculating bias
# parameters
# - donor: fasta file of target
# - reference: fasta file of reference
# - fastq: fastq files
#####################################################

import argparse
import os
import random

import bias
#import config

BAM_TO_SAM="samtools view -h %s"
BEDTOOLS="bedtools"
BOWTIE_PATH="bowtie2"
BWA_PATH="/mnt/work/reference-bias/tools/bwa-0.7.12/bwa"
#BWA_PATH="bwa"
MAUVE_PATH="progressiveMauve"
SAMTOOLS="samtools"

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('--donor', help='donor fasta')
  parser.add_argument('--reference', help='reference fasta')
  parser.add_argument('--job', required=False, help='use to continue previous pipeline')
  parser.add_argument('--start', required=False, help='start from this stage')
  parser.add_argument('--tmpdir', required=False, help='where to write files')
  parser.add_argument('--align', required=False, default='bwa', help='aligner to use')
  parser.add_argument('--donorbam', required=False, help='specify a previously aligned donor bam')
  parser.add_argument('--donorsam', required=False, help='specify a previously aligned donor sam')
  parser.add_argument('fastq', help='fastq files to align')
  args = parser.parse_args()
  # now do each stage...
  calculator = bias.Calculator( BWA_PATH, BOWTIE_PATH, MAUVE_PATH, BAM_TO_SAM, sys.stderr, sys.stdout )
  calculator.calculate( args.donor, args.reference, args.job, args.start, args.tmpdir, args.align, args.donorbam, args.donorsam, args.fastq )
