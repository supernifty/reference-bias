#!/usr/bin/env python

# given an xmfa, remap a bam file coordinates
# usage remap_bam --xmfa xmfa --origin 0 --target 1 srcbamfile > targetsamfile
# writes coverage stats to stderr

import argparse
#import config
import sys

import bias

def remap_bam( xmfa, origin, target, output, new_reference, remap_cigar, output_not_covered, output_target_coverage, bam, out_fh, bam_to_sam ):
  # measure coverage
  print "building coverage map with new_reference {0}...".format( new_reference )
  #print args.origin, dest
  mauve_map = bias.MauveMap( open(xmfa, 'r'), src_strand=origin, target_strand=target, new_reference=new_reference )
  print len(mauve_map.coverage), "positions mapped"

  if output_target_coverage:
    print "generating target bed..."
    target = set()
    for k, v in mauve_map.coverage.iteritems():
      target.add( v )
    print "writing target bed..."
    with open( output_target_coverage, 'w' ) as ofh:
      range_start = None
      for x in xrange(mauve_map.genome_stats['ymin'], mauve_map.genome_stats['ymax']+2): # +2 to make sure range ends
        if x in target:
          if range_start is None:
            range_start = x
        else:
          if range_start is not None:
            ofh.write( '%s\t%i\t%i\n' % ( new_reference, range_start, x ) )
            range_start = None
  
  print "analyzing bam {0}...".format( bam )
  if bam.endswith( '.bam' ):
    sam_fh = bias.BamReaderExternal( bam_to_sam, bam )
  else:
    sam_fh = open( bam, 'r' )
  
  if output_not_covered is not None:
      mauve_map.remap( sam_fh, open( output, 'w' ), remap_cigar=remap_cigar, not_covered_output=open( output_not_covered, 'w' ) )
  else:
      mauve_map.remap( sam_fh, open( output, 'w' ), remap_cigar=remap_cigar )
   
  out_fh.write( "====== Mapping Stats =====\n" )
  for key in mauve_map.genome_stats.keys():
    out_fh.write( "%20s: %i\n" % ( key, mauve_map.genome_stats[key] ) )
  
  out_fh.write( "====== Coverage Stats =====\n" )
  for key in mauve_map.stats.keys():
    out_fh.write( "%20s: %i\n" % ( key, mauve_map.stats[key] ) )
 
if __name__ == '__main__':
  args = parser.parse_args()
  parser = argparse.ArgumentParser(description='Remap a bam file using XMFA')
  parser.add_argument('bam', help='bam file to analyze')
  parser.add_argument('--xmfa', required=True, help='xmfa file')
  parser.add_argument('--origin', type=int, default=1, help='index of fasta that bam was mapped to')
  parser.add_argument('--target', type=int, default=2, help='index of fasta that bam should be mapped to')
  parser.add_argument('--output', required=True, help='sam output file')
  parser.add_argument('--new_reference', required=False, help='new reference name for sam output')
  parser.add_argument('--remap_cigar', required=False, default=False, help='attempt to remap the cigar string')
  parser.add_argument('--output_not_covered', required=False, help='sam output file for aligned reads with no map')
  parser.add_argument('--output_target_coverage', required=False, help='target map')

  remap_bam( args.xmfa, args.origin, args.target, args.output, args.new_reference, args.remap_cigar, args.output_not_covered, args.output_target_coverage, args.bam, sys.stdout, "samtools view -h %s" )
 
