
import datetime
import sys

def log_stderr(msg):
  when = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")
  sys.stderr.write( '%s: %s\n' % ( when, msg ) )

def log_quiet(msg):
  pass

