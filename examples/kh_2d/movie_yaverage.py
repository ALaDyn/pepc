#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import fieldblob as fb

def plot_fieldblob(arg):
  fname, fmin, fmax, ftype = arg
  print "plotting: " + fname
  
  plt.ioff()
  n   = fb.n_of_fieldblob(fname)
  t   = fb.t_of_fieldblob(fname)
  y   = fb.field_of_fieldblob(fname)
  if (ftype == 'density'): y = np.abs(y)
  ymean = np.mean(y, axis = 1)
  offset = fb.offset_of_fieldblob(fname)
  extent = fb.extent_of_fieldblob(fname)
  dx = extent / n

  fig = plt.figure()
  ax  = fig.add_subplot(111)
  
  x = np.linspace(
    offset[0] + 0.5 * dx[0],
    offset[0] + extent[0] - 0.5 * dx[0],
    num = ymean.shape[0]
  )

  ax.plot(x, ymean, 'k-')
  ax.set_xlabel('x / lambda_D')
  ax.set_xlim(offset[0], offset[0] + extent[0])
  ax.set_ylim(fmin, fmax)

  ax.text(0.05, 0.95, 'w_{p,e} t = ' + str(t),
    transform = ax.transAxes)

  fig.savefig(
    os.path.join(os.path.dirname(fname),
      'yaverage_' + os.path.basename(fname) + '.png')
  )

def maxabs_of_fieldblob(fname):
  return np.max(np.abs(np.mean(fb.field_of_fieldblob(fname), axis = 1)))  

def print_usage():
  print "Usage: movie_yaverage.py fieldtype fieldblob1 [fieldblob2 ...]"
  print ""
  print "fieldtype: { density | potential }"

if __name__ == '__main__':
  import multiprocessing as mp
  if (len(sys.argv) < 3):
    print_usage()
    sys.exit(1)
  else:
    ftype  = sys.argv[1]
    fnames = sys.argv[2:]
    pool   = mp.Pool()

    # normalization starting at 0
    if (ftype == 'density'):
      fmin = 0.0
      fmax = np.max(pool.map(maxabs_of_fieldblob, fnames))
    elif (ftype == 'potential'):
      # normalization around 0
      fmax = np.max(pool.map(maxabs_of_fieldblob, fnames))
      fmin = -fmax
    else:
      print_usage()
      sys.exit(1)

    print "fmin: {fmin}".format(fmin = fmin)
    print "fmax: {fmax}".format(fmax = fmax)

    pool.map(plot_fieldblob, ( (fn, fmin, fmax, ftype) for fn in fnames))
