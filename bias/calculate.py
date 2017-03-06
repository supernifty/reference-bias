
import datetime
import os
import random

import bias

class Calculator (object):
  def __init__(self, bwa_path, bowtie_path, mauve_path, bam_to_sam, subread_path, log, out ):
    self.bwa_path = bwa_path
    self.bowtie_path = bowtie_path
    self.mauve_path = mauve_path
    self.bam_to_sam = bam_to_sam
    self.subread_path = subread_path
    self.log_err = log
    self.log_out = out

  def run( self, cmd ):
    '''
      run a system command
    '''
    self.log( cmd )
    os.system( cmd )

  def log( self, msg ):
    when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
    self.log_err.write( '%s: %s\n' % ( when, msg ) )

  def write( self, msg ):
    self.log_out.write( '{0}\n'.format( msg ) )

  def filter_sam( self, out_fn, in_fn, chromosome):
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
 
  def calculate( self, donor, reference, job, stage, tmpdir, align, donorbam, donorsam, fastq, remap_donor, remap_reference, reference_chromosome, donor_chromosome, bed ):
    if job:
      idx = int(job)
    else:
      idx = random.randint(1, 1e6)
    if stage:
      start = int(stage)
    else:
      start = 0
    if tmpdir is None:
      tmpdir = '/tmp'
    self.log( 'Job ID: %i, Starting at stage %i' % (idx, start) )
    # TODO error correction (ec)
  
    # fasta indexes
    stage = 1
    if start <= stage:
      self.index( align, donor )
      self.index( align, reference )
      self.log( 'Stage %i: Indexing completed' % stage )
  
    stage += 1 # start of stage 2
    # alignment (aln)
    if start <= stage:
      if donorsam is None:
        if donor_chromosome is None:
          self.align( align, donor, fastq, '{0}/donor{1}.sam'.format( tmpdir, idx ) )
        else:
          self.align( align, donor, fastq, '{0}/donorfull{1}.sam'.format( tmpdir, idx ) )
          # filter on chromosome
          self.filter_sam("{0}/donor{1}.sam".format( tmpdir, idx ), '{0}/donorfull{1}.sam'.format( tmpdir, idx ), donor_chromosome)
                   
          #self.run( 'samtools view -bhS {0}/donorfull{1}.sam | samtools sort -o {0}/donorfull{1}.bam'.format( tmpdir, idx ) )
          #self.run( 'samtools index {0}/donorfull{1}.bam'.format( tmpdir, idx ) )
          #self.run( "samtools view -h {0}/donorfull{1}.bam {2} > {0}/donor{1}.sam".format( tmpdir, idx, donor_chromosome ))
          #self.run( 'samtools view -bhS {0}/donor{1}.sam | samtools sort -o {0}/donor{1}.bam'.format( tmpdir, idx ) )
          #self.run( 'samtools index {0}/donor{1}.bam'.format( tmpdir, idx ) )
      self.log( 'Stage %i: Donor alignment completed' % stage )
  
    stage += 1 # start of stage 3
    if start <= stage:
      if reference_chromosome is None:
        self.align( align, reference, fastq, '{0}/reference{1}.sam'.format( tmpdir, idx ) )
      else:
        self.align( align, reference, fastq, '{0}/referencefull{1}.sam'.format( tmpdir, idx ) )
        self.filter_sam("{0}/reference{1}.sam".format( tmpdir, idx ), '{0}/referencefull{1}.sam'.format( tmpdir, idx ), donor_chromosome)
        #self.run( 'samtools view -bhS {0}/referencefull{1}.sam | samtools sort -o {0}/referencefull{1}.bam'.format( tmpdir, idx ) )
        #self.run( 'samtools index {0}/referencefull{1}.bam'.format( tmpdir, idx ) )
        #self.run( "samtools view -h {0}/referencefull{1}.bam {2} > {0}/reference{1}.sam".format( tmpdir, idx, reference_chromosome ))
        #self.run( 'samtools view -bhS {0}/reference{1}.sam | samtools sort -o {0}/reference{1}.bam'.format( tmpdir, idx ) )
        #self.run( 'samtools index {0}/reference{1}.bam'.format( tmpdir, idx ) )
      #self.run( '%s mem -t 8 %s %s > %s/reference%i.sam' % ( BWA_PATH, reference, fastq, tmpdir, idx ) )
      self.log( 'Stage %i: Reference alignment completed' % stage )
  
    # genome alignment (mauve)
    if remap_donor is None:
      remap_donor = donor
    if remap_reference is None:
      remap_reference = reference

    stage += 1 # 4
    if start <= stage:
      self.run( '%s --output=%s/mauve%i %s %s' % ( self.mauve_path, tmpdir, idx, remap_donor, remap_reference ) )
      self.log( 'Stage %i: Mauve completed' % stage )

    donor_accession = open( remap_donor, 'r' ).readline().strip().split()[0][1:]
    reference_accession = open( remap_reference, 'r' ).readline().strip().split()[0][1:]

    # realignment
    stage += 1 # 5
    if start <= stage:
      #self.run( 'python remap_bam.py --xmfa %s/mauve%i --origin 2 --target 1 --output_not_covered %s/notcovered%i.sam --output %s/remapped%i.sam %s/reference%i.sam --new_reference \'%s\' > %s/remap_bam%i.stats' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, donor_accession, tmpdir, idx ) )
      #self.log( 'python remap_bam.py --xmfa %s/mauve%i --origin 2 --target 1 --output_not_covered %s/notcovered%i.sam --output %s/remapped%i.sam %s/reference%i.sam --new_reference \'%s\' > %s/remap_bam%i.stats' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, donor_accession, tmpdir, idx ) )
      #xmfa, origin, target, output, new_reference, remap_cigar, output_not_covered, bam, out_fh, bam_to_sam
      bias.remap_bam( xmfa='{0}/mauve{1}'.format( tmpdir, idx ), origin=2, target=1, output_not_covered='{0}/notcovered{1}.sam'.format( tmpdir, idx ), output_target_coverage='{0}/mauve_target{1}.bed'.format( tmpdir, idx ), output='{0}/remapped{1}.sam'.format(tmpdir, idx), new_reference=donor_accession, old_reference=reference_accession, remap_cigar=False, bam='{0}/reference{1}.sam'.format( tmpdir, idx ), out_fh=open( '{0}/remap_bam{1}.stats'.format( tmpdir, idx ), 'w' ), bam_to_sam=self.bam_to_sam )
      self.log( 'Stage %i: Remap completed' % stage )
  
    if donorbam is None:
      target_donorbam = '{0}/donor{1}.bam'.format( tmpdir, idx )
    else:
      target_donorbam = donorbam
  
    # convert to bam
    stage += 1 # 6
    if start <= stage:
      if donorbam is None:
        if donor_accession is None:
          self.run( 'samtools view -bS {0}/donor{1}.sam | samtools sort -o {2}'.format( tmpdir, idx, target_donorbam ) )
        else: # filter on chromosome
          self.run( 'samtools view -bS {0}/donor{1}.sam | samtools sort -o {0}/donorprelim{1}.bam'.format( tmpdir, idx ) )
          self.run( 'samtools index {0}/donorprelim{1}.bam'.format( tmpdir, idx ) )
          self.run( 'samtools view -bh {0}/donorprelim{1}.bam "{2}" > {3}'.format(tmpdir, idx, donor_accession, target_donorbam))
      if reference_accession is None:
        self.run( 'samtools view -bS %s/reference%i.sam | samtools sort -o %s/reference%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
      else:
        self.run( 'samtools view -bS {0}/reference{1}.sam | samtools sort -o {0}/refprelim{1}.bam'.format( tmpdir, idx) )
        self.run( 'samtools index {0}/refprelim{1}.bam'.format( tmpdir, idx) )
        self.run( 'samtools view -bh {0}/refprelim{1}.bam "{2}" > {0}/reference{1}.bam'.format(tmpdir, idx, reference_accession, target_donorbam))
      # fix remapped
      l = []
      if donorsam:
        with open( donorsam, 'r' ) as dfh: # use if already have donor bam
          for line in dfh:
            l.append(line)
            if line.startswith('@PG'):
              break
          #l = (dfh.readline(), dfh.readline())
      else:
        with open( '%s/donor%i.sam' % ( tmpdir, idx ), 'r' ) as dfh:
          for line in dfh:
            l.append(line)
            if line.startswith('@PG'):
              break
          #l = (dfh.readline(), dfh.readline())
  
      with open( '%s/remapped%i.head' % ( tmpdir, idx ), 'w' ) as rfh:
        [ rfh.write(x) for x in l ]
        #rfh.write( l[0] )
        #rfh.write( l[1] )

      #self.run( 'cat %s/remapped%i.head %s/remapped%i.sam | samtools view -bS - > %s/remapped%i.bam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
      self.run( 'samtools view -bS %s/remapped%i.sam | samtools sort -o %s/remapped%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
      l = []
      with open( '%s/reference%i.sam' % ( tmpdir, idx ), 'r' ) as dfh:
        for line in dfh:
          l.append(line)
          if line.startswith('@PG'):
            break
        #l = (dfh.readline(), dfh.readline())
      with open( '%s/notcovered%i.head' % ( tmpdir, idx ), 'w' ) as rfh:
        [ rfh.write(x) for x in l ]
        #rfh.write( l[0] )
        #rfh.write( l[1] )
      self.run( 'cat %s/notcovered%i.head %s/notcovered%i.sam | samtools view -bS - | samtools sort -o %s/notcovered%i.bam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
      #self.run( 'samtools index {0}/notcoveredprelim{1}.bam'.format( tmpdir, idx) )
      #self.run( 'samtools view -bS {0}/notcoveredprelim{1}.bam | samtools sort -o {0}/notcovered{1}.bam' % ( tmpdir, idx ) )
      self.log( 'Stage %i: Convert to bam completed' % stage )

    stage += 1 # 7
    if start <= stage:
      # reads in donor
      self.run( 'samtools flagstat {0} > {1}/donorflag{2}.txt'.format( target_donorbam, tmpdir, idx ) )
      self.run( 'samtools flagstat %s/reference%i.bam > %s/referenceflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'samtools flagstat %s/remapped%i.bam > %s/remappedflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'samtools flagstat %s/notcovered%i.bam > %s/notcoveredflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
     # coverage
  
      #self.run( 'bedtools genomecov -ibam %s/donor%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/donor%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'bedtools genomecov -ibam {0} -bga | awk \'$4<1\' > {1}/donornotcovered{2}.bed'.format( target_donorbam, tmpdir, idx ) )
      # summary of regions that aren't covered on the donor
      self.run( 'awk \'$4<1\' %s/donornotcovered%i.bed | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/donor%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      # summary of regions that aren't covered on the reference
      self.run( 'bedtools genomecov -ibam %s/reference%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/reference%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'bedtools genomecov -ibam %s/remapped%i.bam -bga | awk \'$4<1\' | awk \'{ print $3-$2; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/remapped%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'bedtools genomecov -ibam {0} -d | awk \'$3>0\' | awk \'{{ print $3; }}\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > {1}/donorsum{2}.cov'.format( target_donorbam, tmpdir, idx ) )
      self.run( 'bedtools genomecov -ibam %s/remapped%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/remappedsum%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      self.log( 'Stage %i: Coverage analysis completed' % stage )
  
    stage += 1 # 8
    if start <= stage:
      #self.run( 'python compare_bams.py --compare_position True --subset_detail True --mismatch_detail 1 --xmfa %s/mauve%i --origin 2 --target 1 %s/donor%i.bam %s/remapped%i.bam > %s/compare_bams%i.log' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
      # this was previously run using pypy to save memory
      bias.compare_bams( compare_position=True, mapq=-1, subset_detail=False, mismatch_detail=1, xmfa='{0}/mauve{1}'.format( tmpdir, idx ), origin=2, target=1, bams=( target_donorbam, '{0}/remapped{1}.bam'.format( tmpdir, idx ) ), out_fh=open( '{0}/compare_bams{1}.log'.format(tmpdir, idx), 'w' ), bam_to_sam=self.bam_to_sam )
      #self.run( 'python extract_mismatched_reads.py --min_distance 50 %s/remapped%i.bam < %s/compare_bams%i.log > %s/mismatched%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
      #self.run( 'python extract_mismatched_reads.py --min_distance 1 --max_distance 49 %s/remapped%i.bam < %s/compare_bams%i.log > %s/almost%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
      bias.extract_mismatched_reads( min_distance=50, max_distance=1e9, mappable=False, bam='{0}/remapped{1}.bam'.format( tmpdir, idx ), in_fh=open( '{0}/compare_bams{1}.log'.format( tmpdir, idx ), 'r' ), out_fh=open( '{0}/mismatched{1}.sam'.format( tmpdir, idx ), 'w' ), bam_to_sam=self.bam_to_sam )
      bias.extract_mismatched_reads( min_distance=1, max_distance=49, mappable=False, bam='{0}/remapped{1}.bam'.format( tmpdir, idx ), in_fh=open( '{0}/compare_bams{1}.log'.format( tmpdir, idx ), 'r' ), out_fh=open( '{0}/almost{1}.sam'.format( tmpdir, idx ), 'w' ), bam_to_sam=self.bam_to_sam )
      self.run( 'samtools view -bS %s/mismatched%i.sam > %s/mismatched%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'samtools view -bS %s/almost%i.sam > %s/almost%i.bam' % ( tmpdir, idx, tmpdir, idx ) )
      self.log( 'Stage %i: Mismatch analysis completed' % stage )

    stage += 1 # 9
    if start <= stage:
      self.run( 'samtools flagstat %s/mismatched%i.bam > %s/mismatchedflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'samtools flagstat %s/almost%i.bam > %s/almostflag%i.txt' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'bedtools genomecov -ibam %s/mismatched%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/mismatched%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      self.run( 'bedtools genomecov -ibam %s/notcovered%i.bam -d | awk \'$3>0\' | awk \'{ print $3; }\' | datamash sum 1 min 1 max 1 mean 1 pstdev 1 count 1 > %s/notcovered%i.cov' % ( tmpdir, idx, tmpdir, idx ) )
      if os.stat('{0}/donornotcovered{1}.bed'.format( tmpdir, idx )).st_size == 0:
        with open( "{0}/notcovered_overlap{1}.cov".format( tmpdir, idx ), 'w' ) as fh: 
          fh.write( '0\n' )
      else:
        self.run( "bedtools intersect -a %s/donornotcovered%i.bed -b %s/mauve_target%i.bed | awk '{t+=$3-$2;} END {print t;}' > %s/notcovered_overlap%i.cov" % ( tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
  
      self.log( 'Stage %i: Mismatch coverage completed' % stage )
  
    stage += 1
    if start <= stage:
      # find @SQ SN:tiny2  LN:2310
      reflen = self.find_sequence_len( open( '%s/notcovered%i.head' % ( tmpdir, idx ), 'r' ) )
      remapped_donorlen = self.find_sequence_len_sum( open( '%s/remapped%i.head' % ( tmpdir, idx ), 'r' ) )
      self.write( "===== Stats =====" )
      # reads
      self.write( "-- Reads --" )
      for line in open( '%s/donorflag%i.txt' % ( tmpdir, idx ) ):
        fields = line.strip().split(' ')
        if len(fields) > 4 and fields[4] == 'total':
          self.write( "Donor reads total: %s" % fields[0] )
          drt = int(fields[0])
        if len(fields) > 3 and fields[3] == 'mapped':
          self.write( "Donor reads mapped: %s" % fields[0] )
          drm = int(fields[0])
      self.write( "Donor reads %% mapped: %.1f" % ( 100. * drm / max( 1, drt ) ) )
  
      for line in open( '%s/referenceflag%i.txt' % ( tmpdir, idx ) ):
        fields = line.strip().split()
        if len(fields) > 4 and fields[4] == 'total':
          self.write( "Reference reads total: %s" % fields[0] )
          rrt = int(fields[0])
        if len(fields) > 3 and fields[3] == 'mapped':
          self.write( "Reference reads mapped: %s" % fields[0] )
          rrm = int(fields[0])
      self.write( "Reference reads %% mapped: %.1f" % ( 100. * rrm / max( 1, rrt ) ) )
  
      for line in open( '%s/remappedflag%i.txt' % ( tmpdir, idx ) ):
        fields = line.strip().split()
        if len(fields) > 4 and fields[4] == 'total':
          self.write( "Remapped reads total: %s" % fields[0] )
          mrt = int(fields[0])
        if len(fields) > 3 and fields[3] == 'mapped':
          mrm = int(fields[0])
          self.write( "Remapped reads mapped: %i (%.1f%%)" % (mrm, 100. * mrm / max( 1, rrm ) ) )
      self.write( "Remapped reads %% mapped: %.1f" % ( 100. * mrm / max( 1, rrm ) ) )

      xrm = 0
      xrt = 0
      for line in open( '%s/mismatchedflag%i.txt' % ( tmpdir, idx ) ):
        fields = line.strip().split()
        if len(fields) > 4 and fields[4] == 'total':
          self.write( "Mismatched reads total: %s" % fields[0] )
          xrt = int(fields[0])
        if len(fields) > 3 and fields[3] == 'mapped':
          self.write( "Mismatched reads mapped: %s" % fields[0] )
          xrm = int(fields[0])
      self.write( "Mismatched reads %% mapped: %.1f" % ( 100. * xrm / max(1, xrt) ) )
  
      for line in open( '%s/notcoveredflag%i.txt' % ( tmpdir, idx ) ):
        fields = line.strip().split()
        if len(fields) > 4 and fields[4] == 'total':
          nrt = int(fields[0])
          self.write( "Notcovered reads total: %i (%.1f%%)" % (nrt, 100. * nrt / max( 1, rrt ) ) )
        if len(fields) > 3 and fields[3] == 'mapped':
          nrm = int(fields[0])
          self.write( "Notcovered reads mapped: %i (%.1f%%)" % (nrm, 100. * nrm / max( 1, rrm ) ) )
      self.write( "Notcovered reads %% mapped: %.1f" % ( 100. * nrm / max( 1, nrt ) ) )
  
      for line in open( '%s/almostflag%i.txt' % ( tmpdir, idx ) ):
        fields = line.strip().split()
        if len(fields) > 4 and fields[4] == 'total':
          self.write( "Almost correct reads total: %s" % fields[0] )
          art = int(fields[0])
        if len(fields) > 3 and fields[3] == 'mapped':
          self.write( "Almost correct reads mapped: %s" % fields[0] )
          arm = int(fields[0])
      self.write( "Almost correct reads %% mapped: %.1f" % ( 100. * arm / max( 1, art ) ) )
  
      self.write( "\n-- Correctness --" )
      self.write( "Mapped to correct location: %i (%.1f%%)" % ( mrm - xrm - arm, 100. * ( mrm - xrm - arm ) / max( 1, mrm ) ) )
      self.write( "Mapped correctly or within 50bp: %i (%.1f%%)" % ( mrm - xrm, 100. * ( mrm - xrm ) / max( 1, mrm ) ) )
      self.write( "Mapped incorrectly <50bp: %i (%.1f%%)" % ( arm, 100. * arm / max( 1, mrm ) ) )
      self.write( "Mapped incorrectly >50bp: %i (%.1f%%)" % ( xrm, 100. * xrm / max( 1, mrm ) ) )

      # coverage
      self.write( "\n-- Coverage --" )
      df = open( '%s/donor%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      if len(df) == 0:
        df = (0,0,0,0,0,0)
      donor_not_covered = int(df[0])
      try:
        not_covered_overlap = int( open( '%s/notcovered_overlap%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip() )
      except:
        not_covered_overlap = 0
      self.write( "Donor not covered: %i (%.2f%%)" % ( donor_not_covered, 100. * donor_not_covered / max( 1, remapped_donorlen ) ) )
      self.write( "Donor not covered with mauve target: %i (%.2f%%)" % ( not_covered_overlap, 100. * not_covered_overlap / max( 1, donor_not_covered ) ) )
      self.write( "Donor covered: %i (%.2f%%)" % ( remapped_donorlen - int(df[0]), 100. * (remapped_donorlen - int(df[0]) ) / max( 1, remapped_donorlen ) ) )
      self.write( "Donor gaps: %s" % df[5] )
      self.write( "Donor max gap: %s" % df[2] )
      dfs = open( '%s/donorsum%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      self.write( "Donor mean coverage: %s" % dfs[3] )
      self.write( "Donor max coverage: %s" % dfs[2] )
  
      rf = open( '%s/reference%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      if len(rf) == 0:
        rf = (0,0,0,0,0,0)
      self.write( "Reference not covered: %s (%.1f%%)" % ( rf[0], 100. * int(rf[0]) / max( 1, reflen ) ) )
      self.write( "Reference covered: %i (%.1f%%)" % ( reflen - int(rf[0]), 100. * (reflen - int(rf[0])) / max( 1, reflen ) ) )
      self.write( "Reference gaps: %s" % rf[5] )
      self.write( "Reference max gap: %s" % rf[2] )
  
      mf = open( '%s/remapped%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      if len(mf) == 0:
        mf = (0,0,0,0,0,0)
      self.write( "Remapped not covered: %s (%.1f%%)" % (mf[0], 100. * int(mf[0]) / max( 1, remapped_donorlen ) ) )
      self.write( "Remapped covered: %i (%.1f%%)" % (remapped_donorlen - int(mf[0]), 100. * (remapped_donorlen - int(mf[0]) ) / max( 1, remapped_donorlen ) ) )
      self.write( "Remapped gaps: %s" % mf[5] )
      self.write( "Remapped max gap: %s" % mf[2] )
      mfs = open( '%s/remappedsum%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      self.write( "Remapped mean coverage: %s" % mfs[3] )
      self.write( "Remapped max coverage: %s" % mfs[2] )
  
      self.write( "\n-- Remapped incorrectly > 50bp --" )
      xf = open( '%s/mismatched%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      if len(xf) < 5:
        xf = [ '0', '0', '0', '0', '0', '0' ]
      self.write( "Bases affected by mismatch: %s" % xf[5] )
      self.write( "Max mismatch coverage: %s" % xf[2] )
  
      self.write( "\n-- Off target (outside mappable region) --" )
      nf = open( '%s/notcovered%i.cov' % ( tmpdir, idx ), 'r' ).readline().strip().split()
      if len(nf) < 5:
        nf = [ '0', '0', '0', '0', '0', '0' ]
      self.write( "Off target bases: %s" % nf[5] )
      self.write( "Max coverage of off target: %s" % nf[2] )
      self.write( "\n-- Remapping --" )
      remapping_stats = {}
      for line in open( '%s/remap_bam%i.stats' % ( tmpdir, idx ), 'r' ):
        fields = line.strip().split(':')
        if len(fields) >1:
          remapping_stats[fields[0].strip()] = int(fields[1])
  
      self.write( "Mapped bases: %i (%.1f%%)" % (remapping_stats['count'], 100. * remapping_stats['count'] / max( 1, remapped_donorlen ) ) )
      self.write( "Not mapped bases: %i (%.1f%%)" % (remapped_donorlen - remapping_stats['count'], 100. * ( remapped_donorlen - remapping_stats['count'] )/ remapped_donorlen ) )
      self.write( "Mapped blocks: %s" % remapping_stats['blocks'] )
      self.write( "Covered reads: %i (%.1f)" % ( remapping_stats['reads_covered'], 100. * remapping_stats['reads_covered'] / max( 1, remapping_stats['mapped'] ) ) )
      self.write( "Covered partial reads: %i (%.1f)" % ( remapping_stats['reads_covered'] + remapping_stats['reads_partial'], 100. * ( remapping_stats['reads_covered']  + remapping_stats['reads_partial'] )/ remapping_stats['mapped'] ) )
      self.write( "Not mapped reads: %i (%.1f)" % ( remapping_stats['reads_notcovered'] + remapping_stats['reads_partial'], 100. * ( remapping_stats['reads_notcovered']  + remapping_stats['reads_partial'] ) / max( 1, remapping_stats['mapped'] ) ) )
  
      self.write( "\n-- Summary --" )
      coverage_loss = remapped_donorlen - ( reflen - int(rf[0]) )
      self.write( "Donor not covered by direct alignment: %i (%.2f%%)" % ( donor_not_covered, 100. * donor_not_covered / max( 1, remapped_donorlen ) ) )
      self.write( 'Best case loss from reference coverage: %i / %i: %.1f%%' % ( coverage_loss, remapped_donorlen, 100. * coverage_loss / max( 1, remapped_donorlen ) ) )
      self.write( 'Best case loss from remapping: %i / %i: %.1f%%' % ( remapped_donorlen - remapping_stats['count'], remapped_donorlen, 100. * ( remapped_donorlen - remapping_stats['count'] ) / max( 1, remapped_donorlen ) ) )
      self.write( 'Loss after remap coverage: %i / %i: %.1f%%' % ( int(mf[0]), remapped_donorlen, 100. * int(mf[0]) / max( 1, remapped_donorlen ) ) )
      self.write( 'Loss due to remap: %i / %i: %.1f%%' % ( int(mf[0]) - coverage_loss, remapped_donorlen, 100. * ( int(mf[0]) - coverage_loss ) / max( 1, remapped_donorlen ) ) )
      self.write( 'Potential mismatch impact: %i / %i: %.1f%%' % ( int(xf[5]), remapped_donorlen, 100. * int(xf[5]) / max( 1, remapped_donorlen ) ) )
      self.write( 'Off target: %i / %i: %.1f%%' % ( int(nf[5]), remapped_donorlen, 100. * int(nf[5]) / max( 1, remapped_donorlen ) ) )
      self.write( "Donor not covered with mauve target: %i (%.2f%%)" % ( not_covered_overlap, 100. * not_covered_overlap / max( 1, donor_not_covered ) ) )
      bias_low = 100. * ( int(mf[0]) - int(nf[5]) - not_covered_overlap ) / remapped_donorlen
      bias_mid = 100. * ( int(mf[0]) - not_covered_overlap ) / remapped_donorlen
      bias_high = 100. * ( int(mf[0]) + int(xf[5]) ) / remapped_donorlen
      self.write( 'ESTIMATED BIAS: %.1f -> %.1f -> %.1f' % ( bias_low, bias_mid, bias_high ) )
      self.write( "===== " )
  
      self.log( 'Stage %i: Finished' % stage )
  
    stage += 1
    if start <= stage:
      #self.run( 'rm %s/reference%i.sam %s/mismatched%i.sam %s/remapped%i.sam %s/notcovered%i.sam %s/donor%i.sam' % ( tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx, tmpdir, idx ) )
      #self.run( 'rm {0}/*'.format( tmpdir ) )
      self.log( 'Stage %i: Cleanup finished' % stage )

    return bias_low, bias_mid, bias_high
  
  def index( self, mapper, fasta ):
    if mapper == 'bwa':
      self.run( '%s index %s' % ( self.bwa_path, fasta ) )
    elif mapper == 'bowtie2':
      self.run( '%s-build %s %s-bt2' % ( self.bowtie_path, fasta, fasta ) )
    elif mapper == 'subread':
      self.run( '%s/subread-buildindex -o %s.subidx %s' % ( self.subread_path, fasta, fasta ) )
    else:
      raise "Unsupported aligner"
  
  def align( self, mapper, fasta, fastq, sam ):
    if mapper == 'bwa':
      self.run( '%s mem -t 8 %s %s > %s' % ( self.bwa_path, fasta, fastq, sam ) )
    elif mapper == 'bowtie2':
      self.run( '%s --local -p 16 -x %s-bt2 -U %s --quiet -S %s' % ( self.bowtie_path, fasta, fastq, sam ) ) # -t adds time
    elif mapper == 'subread':
      self.run( '%s/subread-align -t 1 --SAMoutput -i %s.subidx -r %s -o %s' % ( self.subread_path, fasta, fastq, sam ) ) # -t adds time
    else:
      raise "Unsupported aligner"
  
  def find_sequence_len_sum( self, fh ):
    total = 0
    # look for @SQ SN:tiny2  LN:2310
    for line in fh:
      for fragment in line.strip().split():
        if fragment.startswith('LN:'):
          total += int( fragment.split(':',2)[-1] )
    return total  

  def find_sequence_len( self, fh ):
    # look for @SQ SN:tiny2  LN:2310
    for line in fh:
      for fragment in line.strip().split():
        if fragment.startswith('LN:'):
          return int( fragment.split(':',2)[-1] )
    # not found
    return 0
