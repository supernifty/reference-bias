#!/usr/bin/env python
#
# generates a chart illustrating the bias as generated from the calculate_bias command
# also requires a mapping from donor to the label to display on the chart
#
# usage
# =====
# python draw_bias.py genome_map result_files > chart.pdf
# result files are expected to end in accession.fasta
#
# genome_map has the format:
# donor_filename,label,category

import sys

import matplotlib.pyplot as plt
import numpy as np
import pylab
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})

def draw_bias( src, targets, out_fh ):
  #pylab.rcParams['figure.figsize'] = 8, 6
  pylab.rcParams['figure.figsize'] = 6, 6
  fig = plt.figure()
  ax = fig.add_subplot(111)

  color_list = ( '#ff0000', '#603090', '#903030', '#f0a020', '#9090f0', '#f0b0c0', '#b06060', '#f03030', '#30f030', '#f09060', '#607060', '#d0d0d0' )
  path_list = ( '', 'EAEC', 'STEC', 'ExPEC', 'AIEC', 'EPEC', 'ETEC', 'EHEC', 'Commensal', 'Other', 'Shigellosis', 'Unclassified' )

  info = {}
  for line in open(src, 'r'):
    fields = line.strip().split(',')
    if len(fields) > 2:
      accession = fields[0].split('/')[-1][:-6]
      info[accession] = { 'n': fields[1], 'p': fields[2] }
      sys.stderr.write( 'added {0}\n'.format( accession ) )

  records = []
  for target in targets:
    #print target
    fn_fields = target.split('/')[-1].split( '.' )
    accession = '.'.join( [ fn_fields[2], fn_fields[3] ] )
    for line in open( target, 'r' ):
      if line.startswith( 'ESTIMATED' ):
        fields = line.strip().split()
        low = float( fields[2] )
        mid = float( fields[4] )
        high = float( fields[6] )
    #print 'path', info[accession]['p']
    if info[accession]['p'] == 'exclude':
      sys.stderr.write( 'skipped {0}\n'.format( accession ) )
    else:
      record = { 'a': accession, 'l': low, 'm': mid, 'h': high, 'p': info[accession]['p'], 'c': color_list[path_list.index(info[accession]['p'])], 'n': info[accession]['n'] }
      records.append( record )

  #default_color = 'b'
  #midcolor = 'r'
  rangecolor = '#e0e0e0'

  y = 0
  tick_size = 30
  mid_height = 10
  mid_width = 8
  range_height = 4
  y_labels = []
  first = set()
  records = sorted( records, key=lambda k: k['m'] )
  for record in records:
    if record['p'] not in first:
      first.add( record['p'] )
      ax.hlines(y=y, xmin=float(record['l']), xmax=float(record['h']), color=rangecolor, lw=range_height)
      ax.vlines(x=float(record['m']), ymin=y-mid_height, ymax=y+mid_height, color=record['c'], lw=mid_width, label=record['p'])
    else:
      ax.hlines(y=y, xmin=float(record['l']), xmax=float(record['h']), color=rangecolor, lw=range_height)
      ax.vlines(x=float(record['m']), ymin=y-mid_height, ymax=y+mid_height, color=record['c'], lw=mid_width)
    y += tick_size
    y_labels.append( record['n'] )
    sys.stderr.write( 'added record {0}\n'.format( record ) )
  ax.set_ylim(ymin=-tick_size / 2, ymax=y)
  ax.set_xlim(xmin=0)
  ax.tick_params(axis='y', labelsize=6)
  ax.yaxis.tick_right()
  ax.yaxis.set_label_position("right")
  ax.set_yticks(np.arange(0, y, tick_size))
  ax.set_yticklabels(y_labels)
  ax.set_xlabel('Estimated bias (%)' )
  if len(first) > 1:
    leg = ax.legend(loc='lower right', prop={'size':9})
    leg.get_frame().set_alpha(0.8)

  fig.savefig(out_fh, format='pdf', dpi=1000)
  #fig.savefig('%s/%s.png' % (REPORT_DIRECTORY, target), format='png', dpi=300)
  #pylab.rcParams['figure.figsize'] = 6, 4

def genome_sort( gmap ):
  def sorter( x ):
    return gmap[x]
  return sorter

if __name__ == '__main__':
  draw_bias( sys.argv[1], sys.argv[2:], sys.stdout )
