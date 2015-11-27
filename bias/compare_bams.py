
import argparse
import sys

import bias
#import config

def compare_bams( bams, mapq, compare_position, subset_detail, mismatch_detail, xmfa, origin, target, out_fh, bam_to_sam ):
  diff = bias.SamDiff( [ bias.BamReaderExternal( bam_to_sam, sam_file ) for sam_file in bams ], mapq_min=mapq, compare_position=compare_position, subset_detail=subset_detail, mismatch_detail=None if mismatch_detail == -1 else mismatch_detail )

  out_fh.write( "mapq stats\n==========\n" )
  out_fh.write( "i:\tn\tmax\tmin\tmean\tsd\tfilename\n" )
  for idx, stats in enumerate( diff.mapq_stats ):
    out_fh.write( '%i:\t%i\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n' % ( idx, stats['mapped'], stats['max'], stats['min'], stats['mean'], stats['sd'], bams[idx] ) )

  out_fh.write( "\nmapped vs unmapped commonality\n===================\n" )
  for key in sorted( diff.totals.keys() ):
    out_fh.write( ( "{0:0%ib}: {1}\n" % ( len(bams) ) ).format( key, diff.totals[key] ) )

  if compare_position:
    #out_fh.write( "\nmapped vs unmapped commonality including position differences\n===================\n" )
    #for key in sorted( diff.position_totals.keys() ):
    #  out_fh.write( ( "{0:0%ib}: {1}\n" % ( len(bams) ) ).format( key, diff.position_totals[key] ) )
    pass

  if subset_detail:
    out_fh.write( "\nmapq vs position differences\n===================\n" )
    out_fh.write( "i:\tmax\tmin\tmean\tsd\thist\n" )
    for key, value in diff.mapq_subset_stats.items():
      bin_key = ( "{0:0%ib}" % ( len(bams) ) ).format( key )
      out_fh.write( '%s:\t%.2f\t%.2f\t%.2f\t%.2f\t%s\n' % ( bin_key, value['max'], value['min'], value['mean'], value['sd'], value['hist'] ) )

  if mismatch_detail > -1:
    out_fh.write( "\nmismatch details\n===================\n" )
    if xmfa is None:
      out_fh.write( "pos,alt,read_id,wrongness\n" )
    else:
      mauve_map = bias.MauveMap( open(xmfa, 'r'), src_strand=origin, target_strand=target )
      out_fh.write( "pos,alt,read_id,wrongness,nearestmap\n" )
    mismatch_count = 0
    unpaired_count = 0
    for read, value in diff.mismatch_stats.items():
      if 'p' in value and 'a' in value:
        if xmfa is None:
          out_fh.write( '%i,%i,%s,%i\n' % ( value['p'], value['a'], read, value['p'] - value['a'] ) )
        else:
          nearest = mauve_map.find_nearest_target( int(value['a']) )
          out_fh.write( '%i,%i,%s,%i,%i\n' % ( value['p'], value['a'], read, value['p'] - value['a'], nearest ) )
        mismatch_count += 1
      else:
        unpaired_count += 1
        #print "missing values", read, value
    bias.log_stderr( "%i mismatches with incorrect alternatives; %i unpaired reads" % ( mismatch_count, unpaired_count ) )

if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='Compare BAMs')
  parser.add_argument('bams', metavar='bams', nargs='+', help='bam files to analyze')
  parser.add_argument('--mapq', dest='mapq_min', type=int, default=-1, help='only consider alignments with at least this mapq value')
  parser.add_argument('--compare_position', dest='compare_position', type=bool, default=False, help='compare position of alignments')
  parser.add_argument('--subset_detail', dest='subset_detail', type=bool, default=False, help='more detail about subsets')
  parser.add_argument('--mismatch_detail', dest='mismatch_detail', type=int, default=-1, help='more detail about mismatches for this bam index')
  parser.add_argument('--xmfa', dest='xmfa', help='include mapping nearness in mismatch detail')
  parser.add_argument('--origin', type=int, default=1, help='index of fasta that bam was mapped to')
  parser.add_argument('--target', type=int, default=2, help='index of fasta that bam should be mapped to')

  args = parser.parse_args()

  compare_bams( args.bams, args.mapq, args.compare_position, args.subset_detail, args.mismatch_detail, args.xmfa, args.origin, args.target, sys.stdout, "samtools view -h %s" )

