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
import sys

import bias
#import config

BAM_TO_SAM="samtools view -h %s"
BEDTOOLS="bedtools"
BOWTIE_PATH="bowtie2"
BWA_PATH="/mnt/work/reference-bias/tools/bwa-0.7.12/bwa"
#BWA_PATH="bwa"
MAUVE_PATH="progressiveMauve"
SAMTOOLS="samtools"
SUBREAD_PATH="/mnt/work/reference-bias/tools/subread-1.5.0-p3-Linux-x86_64/bin"

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('--donor', help='donor fasta to align to')
  parser.add_argument('--reference', help='reference fasta to align to')
  parser.add_argument('--job', required=False, help='use to continue previous pipeline')
  parser.add_argument('--start', required=False, help='start from this stage')
  parser.add_argument('--tmpdir', required=False, help='where to write files')
  parser.add_argument('--align', required=False, default='bwa', help='aligner to use')
  parser.add_argument('--donorbam', required=False, help='specify a previously aligned donor bam')
  parser.add_argument('--donorsam', required=False, help='specify a previously aligned donor sam')
  parser.add_argument('--remap_donor', required=False, help='donor fasta to remap')
  parser.add_argument('--remap_reference', required=False, help='reference fasta to remap')
  parser.add_argument('fastq', help='fastq files to align')
  args = parser.parse_args()
  # now do each stage...
  calculator = bias.Calculator( BWA_PATH, BOWTIE_PATH, MAUVE_PATH, BAM_TO_SAM, SUBREAD_PATH, sys.stderr, sys.stdout )
  calculator.calculate( args.donor, args.reference, args.job, args.start, args.tmpdir, args.align, args.donorbam, args.donorsam, args.fastq, args.remap_donor, args.remap_reference )
