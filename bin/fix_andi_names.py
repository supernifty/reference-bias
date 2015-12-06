
# change andi name to accession
# fix_names accessions < andi.dist > andi.fixed.dist
# see draw_phylogeny for details

import re
import sys

def shorten( name ):
  short_name = re.sub( '(Escherichia coli|strain|whole genome shotgun sequence|complete sequence|DNA|str\\.|complete genome|chromosome|,|\'|E2348/69 )', '', name.split( ' ', 1 )[1].strip() ).strip()
  short_name = re.sub( '  *', ' ', short_name )
  return short_name

#lines = open('ecoli65.list', 'r').readlines()
lines = open(sys.argv[1], 'r').readlines()
l = sys.stdin.readline()
r = []
out = []
for line in sys.stdin: # andi distances
  fields = line.strip().split(' ')
  found = False
  for x in lines:
    if fields[0] in x:
      #new_name = re.sub( '(Escherichia coli|strain|whole genome shotgun sequence|complete sequence|DNA|str\\.|complete genome|chromosome|,|\'|E2348/69 )', '', x.split( ' ', 1 )[1].strip() ).strip()
      #new_name = re.sub( '  *', ' ', new_name )
      new_name = shorten( x )
      fields[0] = new_name
      r.append( new_name )
      found = True
      break
  if not found:
    sys.stderr.write( 'WARN: {0}\n'.format( fields[0] ) )
  out.append( '%s\n' % ( ','.join(fields) ) )

sys.stdout.write( ',%s\n' % ( ','.join(r) ) ) 
for line in out:
  sys.stdout.write( line )
