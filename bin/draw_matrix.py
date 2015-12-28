#!/usr/bin/env python
#
# generates a chart illustrating the bias as generated from the calculate_bias_matrix command
# also requires a mapping from donor to the label to display on the chart
# usage
# python draw_matrix.py genome_map < matrix.out > chart.pdf
# genome_map has the format:
# donor_filename,label,extra

import sys

import matplotlib.pyplot as plt
import numpy as np
import pylab
from matplotlib import rcParams

rcParams.update({'figure.autolayout': True})
#pylab.rcParams['figure.figsize'] = 8, 6
pylab.rcParams['figure.figsize'] = 12, 8

def genome_sort( gmap ):
  def sorter( x ):
    return gmap[x]
  return sorter

def draw_matrix( in_fh, genome_map, out_fh ):
  gmap = {}
  for line in genome_map:
    x, y, rest = [ z.strip() for z in line.split(',') ]
    gmap[x] = y

  genomes = set()
  data = {}
  first = True
  for line in in_fh:
    sys.stderr.write( line )
    if first: # header
      first = False
      continue
    fields = line.strip().split()
    donor = fields[0]
    reference = fields[1]
    bias = float( fields[4] )
    if donor not in genomes:
      genomes.add( donor )
    if reference not in genomes:
      genomes.add( reference )
    data[ '{0},{1}'.format( donor, reference ) ] = bias

  genome_list = sorted( list( genomes ), key=genome_sort( gmap ) )

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.set_aspect(1)
  ax.set_ylabel('Donor')
  ax.set_xlabel('Reference')

  data_array = np.empty([len( genome_list ), len( genome_list )])
  for xi, xv in enumerate( genome_list ):
    for yi, yv in enumerate( genome_list ):
      data_array[yi][xi] = data[ '{0},{1}'.format( yv, xv ) ]
      ax.annotate(str(data_array[yi][xi]), xy=(xi, yi), horizontalalignment='center', verticalalignment='center')

  res = ax.imshow(data_array, interpolation='nearest', cmap=plt.cm.RdYlGn_r)
    
  mapped_genome_list = [ gmap[x] for x in genome_list ]
  sys.stderr.write( '{0}\n'.format( mapped_genome_list ) )
  ax.set_xticks(np.arange(len(mapped_genome_list)), minor=False)
  ax.set_yticks(np.arange(len(mapped_genome_list)), minor=False)
  ax.set_xticklabels( mapped_genome_list, rotation=45 )
  ax.set_yticklabels( mapped_genome_list )

  cb = fig.colorbar(res, label='% loss')

  fig.savefig( out_fh, format='pdf', dpi=1000)

if __name__ == '__main__':
  draw_matrix( sys.stdin, open( sys.argv[1], 'r' ), sys.stdout )
