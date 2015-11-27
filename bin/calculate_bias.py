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
BWA_PATH="bwa"
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

def align( mapper, fasta, fastq, sam ):
  if mapper == 'bwa':
    run( '%s mem -t 8 %s %s > %s' % ( BWA_PATH, donor, fastq, sam ) )
  elif mapper == 'bowtie2':
    run( '%s --local -p 16 -x %s-bt2 -U %s --quiet -S %s' % ( BOWTIE_PATH, fasta, fastq, sam ) ) # -t adds time

def find_sequence_len( fh ):
  # look for @SQ SN:tiny2  LN:2310
  for line in fh:
    for fragment in line.strip().split():
      if fragment.startswith('LN:'):
        return int( fragment.split(':',2)[-1] )
  # not found
  return 0

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
  if args.job:
    idx = int(args.job)
  else:
    idx = random.randint(1, 1e6)
  if args.start:
    start = int(args.start)
  else:
    start = 0
  if args.tmpdir:
    tmpdir = args.tmpdir
  else:
    tmpdir = '/tmp'
  bias.log_stderr( 'Job ID: %i, Starting at stage %i' % (idx, start) )
  # TODO error correction (ec)

  # fasta indexes
  stage = 1
  if start <= stage:
    index( args.align, args.donor )
    index( args.align, args.reference )
    bias.log_stderr( 'Stage %i: Indexing completed' % stage )

  stage += 1 # 2
  # alignment (aln)
  if start <= stage:
    if not args.donorsam:
      align( args.align, args.donor, args.fastq, '{0}/donor{1}.sam'.format( tmpdir, idx ) )
    bias.log_stderr( 'Stage %i: Donor alignment completed' % stage )

  stage += 1 # 3
  if start <= stage:
    align( args.align, args.reference, args.fastq, '{0}/reference{1}.sam'.format( tmpdir, idx ) )
    #run( '%s mem -t 8 %s %s > %s/reference%i.sam' % ( BWA_PATH, args.reference, args.fastq, tmpdir, idx ) )
    bias.log_stderr( 'Stage %i: Reference alignment completed' % stage )

  # genome alignment (mauve)
  stage += 1 # 4
  if start <= stage:
    run( '%s --output=%s/mauve%i %s %s' % ( MAUVE_PATH, tmpdir, idx, args.donor, args.reference ) )
    bias.log_stderr( 'Stage %i: Mauve completed' % stage )

  # realignment
  stage += 1 # 5
  if start <= stage:
    donor_accession = open( args.donor, 'r' ).readline().strip().split()[0][1:]
    #run( 'python remap_bam.py --xmfa %s/mauve%i --origin 2 --target 1 --output_not_covered %s/notcovered%i.sam --output %s/remapped%i.sam %s/reference%i.sam --new_reference \'%s\' > %s/remap_bam%i.stats' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, donor_accession, tmpdir, idx ) )
    bias.log_stderr( 'python remap_bam.py --xmfa %s/mauve%i --origin 2 --target 1 --output_not_covered %s/notcovered%i.sam --output %s/remapped%i.sam %s/reference%i.sam --new_reference \'%s\' > %s/remap_bam%i.stats' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, donor_accession, tmpdir, idx ) )
    #xmfa, origin, target, output, new_reference, remap_cigar, output_not_covered, bam, out_fh, bam_to_sam
    bias.remap_bam( xmfa='{0}/mauve{1}'.format( tmpdir, idx ), origin=2, target=1, output_not_covered='{0}/notcovered{1}.sam'.format( tmpdir, idx ), output_target_coverage='{0}/mauve_target{1}.bed'.format( tmpdir, idx ), output='{0}/remapped{1}.sam'.format(tmpdir, idx), new_reference=donor_accession, remap_cigar=False, bam='{0}/reference{1}.sam'.format( tmpdir, idx ), out_fh=open( '{0}/remap_bam{1}.stats'.format( tmpdir, idx ), 'w' ), bam_to_sam=BAM_TO_SAM )
    bias.log_stderr( 'Stage %i: Remap completed' % stage )

  # convert to bam
  stage += 1 # 6
  if start <= stage:
    if args.donorbam:
      run( 'ln -s {0} {1}/donor{2}.bam' % ( args.donorbam, tmpdir, idx ) ) # use if already have donor bam
    else:
      run( 'samtools view -bS %s/donor%i.sam > %s/donor%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools view -bS %s/reference%i.sam > %s/reference%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    # fix remapped
    if args.donorsam:
      with open( '%s/donor.sam' % ( tmpdir, idx ), 'r' ) as dfh: # use if already have donor bam
        l = (dfh.readline(), dfh.readline())
    else:
      with open( '%s/donor%i.sam' % ( tmpdir, idx ), 'r' ) as dfh:
        l = (dfh.readline(), dfh.readline())

    with open( '%s/remapped%i.head' % ( tmpdir, idx ), 'w' ) as rfh:
      rfh.write( l[0] )
      rfh.write( l[1] )
    #run( 'cat %s/remapped%i.head %s/remapped%i.sam | samtools view -bS - > %s/remapped%i.bam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    run( 'samtools view -bS %s/remapped%i.sam > %s/remapped%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    with open( '%s/reference%i.sam' % ( tmpdir, idx ), 'r' ) as dfh:
      l = (dfh.readline(), dfh.readline())
    with open( '%s/notcovered%i.head' % ( tmpdir, idx ), 'w' ) as rfh:
      rfh.write( l[0] )
      rfh.write( l[1] )
    run( 'cat %s/notcovered%i.head %s/notcovered%i.sam | samtools view -bS - > %s/notcovered%i.bam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    bias.log_stderr( 'Stage %i: Convert to bam completed' % stage )

  stage += 1 # 7
  if start <= stage:
    # reads in donor
    run( 'samtools flagstat %s/donor%i.bam > %s/donorflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/reference%i.bam > %s/referenceflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/remapped%i.bam > %s/remappedflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/notcovered%i.bam > %s/notcoveredflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
   # coverage
    
    #run( 'bedtools genomecov -ibam %s/donor%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/donor%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/donor%i.bam -bga | awk \'$4<1\' > %s/donornotcovered%i.bed' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'awk \'$4<1\' %s/donornotcovered%i.bed | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/donor%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/reference%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/reference%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/remapped%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/remapped%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/donor%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/donorsum%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/remapped%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/remappedsum%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    bias.log_stderr( 'Stage %i: Coverage analysis completed' % stage )

  stage += 1 # 8
  if start <= stage:
    #run( 'python compare_bams.py --compare_position True --subset_detail True --mismatch_detail 1 --xmfa %s/mauve%i --origin 2 --target 1 %s/donor%i.bam %s/remapped%i.bam > %s/compare_bams%i.log' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    # this was previously run using pypy to save memory
    bias.compare_bams( compare_position=True, mapq=-1, subset_detail=False, mismatch_detail=1, xmfa='{0}/mauve{1}'.format( tmpdir, idx ), origin=2, target=1, bams=( '{0}/donor{1}.bam'.format( tmpdir, idx ), '{0}/remapped{1}.bam'.format( tmpdir, idx ) ), out_fh=open( '{0}/compare_bams{1}.log'.format(tmpdir, idx), 'w' ), bam_to_sam=BAM_TO_SAM )
    #run( 'python extract_mismatched_reads.py --min_distance 50 %s/remapped%i.bam < %s/compare_bams%i.log > %s/mismatched%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    #run( 'python extract_mismatched_reads.py --min_distance 1 --max_distance 49 %s/remapped%i.bam < %s/compare_bams%i.log > %s/almost%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    bias.extract_mismatched_reads( min_distance=50, max_distance=1e9, mappable=False, bam='{0}/remapped{1}.bam'.format( tmpdir, idx ), in_fh=open( '{0}/compare_bams{1}.log'.format( tmpdir, idx ), 'r' ), out_fh=open( '{0}/mismatched{1}.sam'.format( tmpdir, idx ), 'w' ), bam_to_sam=BAM_TO_SAM )
    bias.extract_mismatched_reads( min_distance=1, max_distance=49, mappable=False, bam='{0}/remapped{1}.bam'.format( tmpdir, idx ), in_fh=open( '{0}/compare_bams{1}.log'.format( tmpdir, idx ), 'r' ), out_fh=open( '{0}/almost{1}.sam'.format( tmpdir, idx ), 'w' ), bam_to_sam=BAM_TO_SAM )
    run( 'samtools view -bS %s/mismatched%i.sam > %s/mismatched%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools view -bS %s/almost%i.sam > %s/almost%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
    bias.log_stderr( 'Stage %i: Mismatch analysis completed' % stage )

  stage += 1 # 9
  if start <= stage:
    run( 'samtools flagstat %s/mismatched%i.bam > %s/mismatchedflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'samtools flagstat %s/almost%i.bam > %s/almostflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/mismatched%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/mismatched%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( 'bedtools genomecov -ibam %s/notcovered%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/notcovered%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
    run( "bedtools intersect -a %s/donornotcovered%i.bed -b %s/mauve_target%i.bed | awk '{t+=$3-$2;} END {print t;}' > %s/notcovered_overlap%i.cov" % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    
    bias.log_stderr( 'Stage %i: Mismatch coverage completed' % stage )

  stage += 1
  if start <= stage:
    # find @SQ SN:tiny2  LN:2310
    reflen = find_sequence_len( open( '%s/notcovered%i.head' % ( tmpdir, idx ), 'r' ) ) 
    donorlen = find_sequence_len( open( '%s/remapped%i.head' % ( tmpdir, idx ), 'r' ) )
    print "===== Stats ====="
    # reads
    print "-- Reads --"
    for line in open( '%s/donorflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split(' ')
      if len(fields) > 4 and fields[4] == 'total':
        print "Donor reads total: %s" % fields[0]
        drt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Donor reads mapped: %s" % fields[0]
        drm = int(fields[0])
    print "Donor reads %% mapped: %.1f" % ( 100. * drm / max( 1, drt ) )

    for line in open( '%s/referenceflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Reference reads total: %s" % fields[0]
        rrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Reference reads mapped: %s" % fields[0]
        rrm = int(fields[0])
    print "Reference reads %% mapped: %.1f" % ( 100. * rrm / max( 1, rrt ) )

    for line in open( '%s/remappedflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Remapped reads total: %s" % fields[0]
        mrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        mrm = int(fields[0])
        print "Remapped reads mapped: %i (%.1f%%)" % (mrm, 100. * mrm / max( 1, rrm ) )
    print "Remapped reads %% mapped: %.1f" % ( 100. * mrm / max( 1, rrm ) )

    xrm = 0
    xrt = 0
    for line in open( '%s/mismatchedflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Mismatched reads total: %s" % fields[0]
        xrt = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Mismatched reads mapped: %s" % fields[0]
        xrm = int(fields[0])
    print "Mismatched reads %% mapped: %.1f" % ( 100. * xrm / max(1, xrt) )

    for line in open( '%s/notcoveredflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        nrt = int(fields[0])
        print "Notcovered reads total: %i (%.1f%%)" % (nrt, 100. * nrt / max( 1, rrt ) )
      if len(fields) > 3 and fields[3] == 'mapped':
        nrm = int(fields[0])
        print "Notcovered reads mapped: %i (%.1f%%)" % (nrm, 100. * nrm / max( 1, rrm ) )
    print "Notcovered reads %% mapped: %.1f" % ( 100. * nrm / max( 1, nrt ) )

    for line in open( '%s/almostflag%i.txt' % ( tmpdir, idx ) ):
      fields = line.strip().split()
      if len(fields) > 4 and fields[4] == 'total':
        print "Almost correct reads total: %s" % fields[0]
        art = int(fields[0])
      if len(fields) > 3 and fields[3] == 'mapped':
        print "Almost correct reads mapped: %s" % fields[0]
        arm = int(fields[0])
    print "Almost correct reads %% mapped: %.1f" % ( 100. * arm / max( 1, art ) )

    print "\n-- Correctness --"
    print "Mapped to correct location: %i (%.1f%%)" % ( mrm - xrm - arm, 100. * ( mrm - xrm - arm ) / max( 1, mrm ) )
    print "Mapped correctly or within 50bp: %i (%.1f%%)" % ( mrm - xrm, 100. * ( mrm - xrm ) / max( 1, mrm ) )
    print "Mapped incorrectly <50bp: %i (%.1f%%)" % ( arm, 100. * arm / max( 1, mrm ) )
    print "Mapped incorrectly >50bp: %i (%.1f%%)" % ( xrm, 100. * xrm / max( 1, mrm ) )

    # coverage
    print "\n-- Coverage --"
    df = open( '%s/donor%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    if len(df) == 0:
      df = (0,0,0,0,0,0)
    donor_not_covered = int(df[0])
    try:
      not_covered_overlap = int( open( '%s/notcovered_overlap%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip() )
    except:
      not_covered_overlap = 0
    print "Donor not covered: %i (%.2f%%)" % ( donor_not_covered, 100. * donor_not_covered / max( 1, donorlen ) )
    print "Donor not covered with mauve target: %i (%.2f%%)" % ( not_covered_overlap, 100. * not_covered_overlap / max( 1, donor_not_covered ) )
    print "Donor covered: %i (%.2f%%)" % ( donorlen - int(df[0]), 100. * (donorlen - int(df[0]) ) / max( 1, donorlen ) )
    print "Donor gaps: %s" % df[5]
    print "Donor max gap: %s" % df[2]
    dfs = open( '%s/donorsum%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Donor mean coverage: %s" % dfs[3]
    print "Donor max coverage: %s" % dfs[2]

    rf = open( '%s/reference%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Reference not covered: %s (%.1f%%)" % ( rf[0], 100. * int(rf[0]) / max( 1, reflen ) )
    print "Reference covered: %i (%.1f%%)" % ( reflen - int(rf[0]), 100. * (reflen - int(rf[0])) / max( 1, reflen ) )
    print "Reference gaps: %s" % rf[5]
    print "Reference max gap: %s" % rf[2]

    mf = open( '%s/remapped%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Remapped not covered: %s (%.1f%%)" % (mf[0], 100. * int(mf[0]) / max( 1, donorlen ) )
    print "Remapped covered: %i (%.1f%%)" % (donorlen - int(mf[0]), 100. * (donorlen - int(mf[0]) ) / max( 1, donorlen ) )
    print "Remapped gaps: %s" % mf[5]
    print "Remapped max gap: %s" % mf[2]
    mfs = open( '%s/remappedsum%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    print "Remapped mean coverage: %s" % mfs[3]
    print "Remapped max coverage: %s" % mfs[2]

    print "\n-- Remapped incorrectly > 50bp --"
    xf = open( '%s/mismatched%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    if len(xf) < 5:
      xf = [ '0', '0', '0', '0', '0', '0' ]
    print "Bases affected by mismatch: %s" % xf[5]
    print "Max mismatch coverage: %s" % xf[2]

    print "\n-- Off target (outside mappable region) --"
    nf = open( '%s/notcovered%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
    if len(nf) < 5:
      nf = [ '0', '0', '0', '0', '0', '0' ]
    print "Off target bases: %s" % nf[5]
    print "Max coverage of off target: %s" % nf[2]

    print "\n-- Remapping --"
    remapping_stats = {}
    for line in open( '%s/remap_bam%i.stats' % ( tmpdir, idx ), 'r' ):
      fields = line.strip().split(':')
      if len(fields) >1:
        remapping_stats[fields[0].strip()] = int(fields[1])

    print "Mapped bases: %i (%.1f%%)" % (remapping_stats['count'], 100. * remapping_stats['count'] / max( 1, donorlen ) )
    print "Not mapped bases: %i (%.1f%%)" % (donorlen - remapping_stats['count'], 100. * ( donorlen - remapping_stats['count'] )/ donorlen )
    print "Mapped blocks: %s" % remapping_stats['blocks']
    print "Covered reads: %i (%.1f)" % ( remapping_stats['reads_covered'], 100. * remapping_stats['reads_covered'] / max( 1, remapping_stats['mapped'] ) )
    print "Covered partial reads: %i (%.1f)" % ( remapping_stats['reads_covered'] + remapping_stats['reads_partial'], 100. * ( remapping_stats['reads_covered']  + remapping_stats['reads_partial'] )/ remapping_stats['mapped'] )
    print "Not mapped reads: %i (%.1f)" % ( remapping_stats['reads_notcovered'] + remapping_stats['reads_partial'], 100. * ( remapping_stats['reads_notcovered']  + remapping_stats['reads_partial'] ) / max( 1, remapping_stats['mapped'] ) )

    print "\n-- Summary --"
    coverage_loss = donorlen - ( reflen - int(rf[0]) )
    print "Donor not covered by direct alignment: %i (%.2f%%)" % ( donor_not_covered, 100. * donor_not_covered / max( 1, donorlen ) )
    print 'Best case loss from reference coverage: %i / %i: %.1f%%' % ( coverage_loss, donorlen, 100. * coverage_loss / max( 1, donorlen ) )
    print 'Best case loss from remapping: %i / %i: %.1f%%' % ( donorlen - remapping_stats['count'], donorlen, 100. * ( donorlen - remapping_stats['count'] ) / max( 1, donorlen ) )
    print 'Loss after remap coverage: %i / %i: %.1f%%' % ( int(mf[0]), donorlen, 100. * int(mf[0]) / max( 1, donorlen ) )
    print 'Loss due to remap: %i / %i: %.1f%%' % ( int(mf[0]) - coverage_loss, donorlen, 100. * ( int(mf[0]) - coverage_loss ) / max( 1, donorlen ) )
    print 'Potential mismatch impact: %i / %i: %.1f%%' % ( int(xf[5]), donorlen, 100. * int(xf[5]) / max( 1, donorlen ) )
    print 'Off target: %i / %i: %.1f%%' % ( int(nf[5]), donorlen, 100. * int(nf[5]) / max( 1, donorlen ) )
    print "Donor not covered with mauve target: %i (%.2f%%)" % ( not_covered_overlap, 100. * not_covered_overlap / max( 1, donor_not_covered ) )
    print 'ESTIMATED BIAS: %.1f -> %.1f -> %.1f' % ( 100. * ( int(mf[0]) - int(nf[5]) - not_covered_overlap ) / donorlen , 100. * ( int(mf[0]) - not_covered_overlap ) / donorlen, 100. * ( int(mf[0]) + int(xf[5]) ) / donorlen )
    print "===== "

    bias.log_stderr( 'Stage %i: Finished' % stage )

  stage += 1
  if start <= stage:
    #run( 'rm %s/reference%i.sam %s/mismatched%i.sam %s/remapped%i.sam %s/notcovered%i.sam %s/donor%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
    bias.log_stderr( 'Stage %i: Cleanup finished' % stage )
  

