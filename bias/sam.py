
import bias
import collections
import numpy
import subprocess

def bucket( values, buckets ):
  '''
    puts values into buckets
    @values: list of values
    @buckets: list of bucket separators with each being the minimum value that can go in that bucket
  '''
  result = [ 0 ] * len(buckets)
  current_bucket = 0
  for value in sorted( values ):
    while current_bucket < len(buckets) - 1 and value >= buckets[current_bucket + 1]:
      current_bucket += 1
    result[current_bucket] += 1
  return result

class SamDiff(object):
  def __init__( self, sam_fhs, mapq_min=-1, compare_position=False, subset_detail=False, mismatch_detail=None, log=bias.log_stderr, mapq_stats=False ):
    self.log = log
    self.use_mapq_stats = mapq_stats
    self.mapq_min = mapq_min
    self.compare_position = compare_position
    self.subset_detail = subset_detail
    self.mismatch_detail = None if mismatch_detail is None else 2 ** mismatch_detail
    self.mismatch_stats = {}
    self._stats = collections.defaultdict(int) # all sam entries with a binary collection showing which sam mapped it
    self._position_stats = collections.defaultdict(int) # all sam entries and positions, containing a binary collection of mapped sams
    self.mapq_stats = [] # overall mapq of a sample
    self.mapq_subsets = collections.defaultdict(list)  # mapq and the count of subsets that it is assigned to i.e. { 15: { '00': 12, '01': 14 ... etc } }

    for idx, reader in enumerate(sam_fhs): # process each sam file
      log( 'processing file %i...' % idx )
      bit_pos = 2 ** idx
      self.mapq = []
      for pos, line in enumerate(reader):
        self.parse_line( pos, bit_pos, line.strip() )
        if ( pos + 1 ) % 100000 == 0:
          log( 'processed %i lines...' % (pos+1) )
      if mapq_stats:
        if len(self.mapq) > 0:
          self.mapq_stats.append( { 'mapped': len(self.mapq), 'max': max(self.mapq), 'min': min(self.mapq), 'mean': numpy.mean( self.mapq ), 'sd': numpy.std( self.mapq ) } )
        else:
          self.mapq_stats.append( { 'mapped': 0, 'max': 0, 'min': 0, 'mean': 0, 'sd': 0 } )
      log( 'processed file %i (bit_pos %i): read %i lines' % ( idx, bit_pos, pos+1 ) )

    log( 'analyzing...' )
    self.totals = collections.defaultdict(int)
    for key, value in self._stats.iteritems(): # readname, distribution
      self.totals[value] += 1
    # free up memory
    del self._stats

    if self.compare_position:
      log( 'analyzing positions...' )
      #self.position_totals = collections.defaultdict(int)
      self.mapq_totals = collections.defaultdict(list)
      for key, value in self._position_stats.iteritems(): # readname-pos, distribution
        #self.position_totals[value] += 1
        if mapq_stats and self.subset_detail:
          self.mapq_totals[value].extend( self.mapq_subsets[key] )
        if self.mismatch_detail is not None and self.mismatch_detail == value: # mismatch_detail exactly
          name, pos = key.split( '-' )
          self.mismatch_stats[name] = { 'p': int(pos) }
          #log( 'found %s at %i' % (name, int(pos)) )
      log( 'analyzing mismatches...' )
      if self.mismatch_detail is not None: # 2nd pass for mismatch updates
        for key, value in self._position_stats.iteritems(): # readname-pos, distribution
          if value != 0 and value & self.mismatch_detail == 0: # mismatch not present
            name, pos = key.split( '-' )
            if name in self.mismatch_stats:
              self.mismatch_stats[name]['a'] = int(pos)
              #log( 'found alt for %s at %i' % (name, int(pos)) )
            else:
              pass
              #log( 'unmatched read %s' % name )
      # free up memory
      del self._position_stats

      # now generate stats
      if mapq_stats:
        log( 'analyzing mapq distribution...' )
        self.mapq_subset_stats = {}
        buckets = (0, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60)
        for key, value in self.mapq_totals.iteritems():
          if len(value) > 0:
            min_val = min(value)
            max_val = max(value)
            self.mapq_subset_stats[key] = { 'max': max_val, 'min': min_val, 'mean': numpy.mean(value), 'sd': numpy.std(value), 'hist': bucket(value, buckets)  }

  def parse_line( self, pos, bit_pos, line ):
    fields = line.split()
    if len(fields) < 5 or line.startswith('@'):
      pass #self.log( 'WARN: %i: unexpected format: %s' % ( pos, line.strip() ) )
    else:
      flag = int( fields[1] )
      mapq = int( fields[4] )
      if flag & 0x04 != 0 or mapq < self.mapq_min: # unmapped
        self._stats[fields[0]] |= 0 # unmapped read - don't change; set to 0 if not already present
      else:
        self.mapq.append( mapq )
        if flag & 0x02 != 0: # mapped read
          self._stats[fields[0]] |= bit_pos
        else: # unknown mapping (assume mapped)
          self._stats[fields[0]] |= bit_pos
      if self.compare_position:
        position = int( fields[3] )
        ident = '%s-%i' % ( fields[0], position )
        if flag & 0x04 != 0 or mapq < self.mapq_min: # unmapped
          self._position_stats[ident] |= 0 # unmapped read - don't change; set to 0 if not already present
        else:
          # position tracking
          if flag & 0x02 != 0: # mapped read
            self._position_stats[ident] |= bit_pos
          else: # unknown mapping (assume mapped)
            self._position_stats[ident] |= bit_pos
          if self.subset_detail:
            self.mapq_subsets[ident].append( mapq )

class SamFilter(object):
  '''
    filter sam reads based on tag names
  '''
  def __init__( self, sam_fh, target_fh, allowed_tags, exclude=False, log=bias.log_stderr ):
    log( 'SamFilter: starting...' )
    self.allowed_tags = allowed_tags
    self.exclude = exclude
    written = 0
    pos = 0
    for pos, line in enumerate(sam_fh):
      if self.parse_line( line.strip() ):
        target_fh.write( '%s\n' % line.strip() )
        written += 1
      if ( pos + 1 ) % 100000 == 0:
        log( 'read %i lines, wrote %i lines...' % (pos+1, written) )
    log( 'processed: read %i lines, wrote %i lines' % ( pos+1, written ) )

  def parse_line( self, line ):
    fields = line.split()
    if len(fields) < 5 or line.startswith('@'):
      return True  
    else:
      if fields[0] in self.allowed_tags:
        return not self.exclude
      else:
        return self.exclude
  
class BamReaderExternal(object):
  '''
    sam handle interface using an external BamToSam
  '''
  def __init__( self, cmd, sam_file ):
    p = subprocess.Popen(cmd % sam_file, shell=True, bufsize=0, stdin=subprocess.PIPE, stdout=subprocess.PIPE) #, close_fds=True)
    self.stdout = p.stdout

  def __iter__(self):
    return self.stdout.__iter__()
  
