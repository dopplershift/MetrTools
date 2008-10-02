import numpy as N

def rms(diffs):
  return N.sqrt(N.average(diffs**2))
