
# to build the phylogeny:
# - andi *.fasta > ecoli.distances
# - python bin/fix_andi_names.py accession_list < ../results/ecoli.distances > ../results/ecoli.distances.fixed # turns andi prefixes into genome names
# - python bin/draw_phylogeny.py ../results/ecoli_matrix_2.map [full.path.to]/ecoli.distances.fixed ../results/EDL933/calculate_bias.EDL933.* > maketree.R
# - run maketree.R to plot the phylogeny
# - also generated is gradient.pdf which can be used as the colour key
#
# usage
# draw_phylogeny genome_map distance_matrix resultfiles > Rcode
#

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import matplotlib.cm
#import matplotlib.colors as mcolors

def parse( src, distances, results ):
  info = {}
  for line in open(src, 'r'):
    fields = line.strip().split(',')
    if len(fields) > 2:
      accession = fields[0].split('/')[-1][:-6]
      info[accession] = { 'n': fields[1], 'p': fields[2] }
      #print 'added', accession, 'with', fields[1]

  records = []
  for target in results:
    #print target
    fn_fields = target.split('/')[-1].split( '.' )
    accession = '.'.join( [ fn_fields[2], fn_fields[3] ] )
    for line in open( target, 'r' ):
      if line.startswith( 'ESTIMATED' ):
        fields = line.strip().split()
        low = float( fields[2] )
        mid = float( fields[4] )
        high = float( fields[6] )
    #print 'name', info[accession]['n'], 'given', accession
    record = { 'a': accession, 'l': low, 'm': mid, 'h': high, 'p': info[accession]['p'], 'n': info[accession]['n'] }
    records.append( record ) 

  return info, records

def make_colour_code( info, distances, records, out_fh ):

  max_bias = 30.
  for record in records:
    max_bias = max( max_bias, record['m'] )
  max_bias = ( max_bias + 5 ) - ( max_bias % 5 ) # round to 5
  out_fh.write( '## max bias is {0}\n'.format( max_bias ) )
  #conv = matplotlib.cm.RdYlGn
  #conv = matplotlib.cm.winter
  conv = matplotlib.colors.LinearSegmentedColormap.from_list(name='gor', colors =['green', 'orange', 'red'], N=16)
  
  # index of data
  #mid_i = fs.index('mid')
  #rnm_i = fs.index('Reference not covered')
  #rl_i = fs.index('Reference Length')

  out_fh.write( 'labelCol <- function(x) {\n' )
  out_fh.write( 'if (is.leaf(x)) {\n' )
  out_fh.write( '## fetch label\n' )
  out_fh.write( 'label <- attr(x, "label") \n' )
  out_fh.write( 'attr(x, "nodePar") <- list(cex=0.5, lab.cex = 0.6, pch=15, bg="#ff0000")\n' )
  out_fh.write( '## set label color to red for A and B, to blue otherwise\n' )

  for record in records:
    name = record['n']
    # data value
    #bias = float(fs[mid_i])
    bias = record['m']
    #new_name = re.sub( '(Escherichia coli|strain|whole genome shotgun sequence|complete sequence|DNA|str\\.|complete genome|chromosome|,|\'|E2348/69 )', '', name ).strip()
    #new_name = re.sub( '  *', ' ', new_name )
    color = matplotlib.colors.rgb2hex( conv( bias / max_bias ) )
    out_fh.write( 'if (label == "%s") { attr(x, "nodePar") <- list(cex=1.0, lab.cex = 0.6, pch = 19, lab.col="%s", bg="%s", col="%s") }\n' % ( name, color, color, color ) )
  
  out_fh.write( '}\n' )
  out_fh.write( 'return(x)\n' )
  out_fh.write( '}\n' )
  # write out clustering code
  out_fh.write( 'distances = read.csv(file = \'{0}\', header = T, row.names= 1)\n'.format( distances ) )
  out_fh.write( 'hc <- hclust(as.dist(distances))\n' )
  out_fh.write( 'd <- dendrapply(as.dendrogram(hc, hang=0.1), labelCol)\n' )
  out_fh.write( 'plot(d, horiz=F)\n' )

  # draw scale
  
  fig = plt.figure()
  ax = fig.add_subplot(111)
  
  a = np.linspace(0, 1, 256).reshape(1,-1)
  a = np.vstack((a,a))
  ax.imshow(a, aspect='auto', cmap=conv)
  
  fig.savefig('gradient.pdf', format='pdf', dpi=1000)

if __name__ == '__main__':
  info, records = parse( sys.argv[1], sys.argv[2], sys.argv[3:] )
  make_colour_code( info, sys.argv[2], records, sys.stdout )
