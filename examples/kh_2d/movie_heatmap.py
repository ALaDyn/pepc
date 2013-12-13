#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import fieldblob as fb

def plot_fieldblob(arg):
  fname, fmin, fmax, ftype, cmname = arg
  print "plotting: " + fname
  dpi = 80.
  
  n   = fb.n_of_fieldblob(fname)
  t   = fb.t_of_fieldblob(fname)
  y   = fb.field_of_fieldblob(fname)
  if (ftype == 'density'): y = np.abs(y)
  fig = plt.figure(figsize = (n[0] / dpi, n[1] / dpi))
  cm  = plt.get_cmap(cmname)

  ax  = fig.add_axes([0,0,1,1], frameon = False)
  
  offset = fb.offset_of_fieldblob(fname)
  extent = fb.extent_of_fieldblob(fname)
  ax.imshow(
    y.transpose(),
    interpolation = 'nearest',
    cmap = cm,
    norm = Normalize(vmin = fmin, vmax = fmax),
    extent = (offset[0], extent[0] + offset[0], offset[1], extent[1] + offset[1]),
    origin = 'bottom'
  )

  ax.set_xticks([])
  ax.set_yticks([])

  ax.text(0.05, 0.95, 'w_{{p,e}} t = {0:6.3f}'.format(t), transform = ax.transAxes),

  fig.savefig(
    os.path.join(os.path.dirname(fname),
      'heatmap_' + os.path.basename(fname) + '.png'),
    dpi = dpi
  )

def maxabs_of_fieldblob(fname):
  return np.max(np.abs(fb.field_of_fieldblob(fname)))

def print_usage():
  print "Usage: movie_heatmap.py fieldtype colormap fieldblob1 [fieldblob2 ...]"
  print ""
  print "fieldtype: { density | potential }"
  print "colormap : any color map known to matplotlib, try 'Blues'/'Reds' or 'RdBu_r'"

if __name__ == '__main__':
  import multiprocessing as mp
  if (len(sys.argv) < 4):
    print_usage()
    sys.exit(1)
  else:
    ftype  = sys.argv[1]
    cmname = sys.argv[2]
    fnames = sys.argv[3:]
    pool   = mp.Pool()

    if (ftype == 'density'): # normalize starting at 0
      fmax = np.max(pool.map(maxabs_of_fieldblob, fnames))
      fmin = 0.0
    elif (ftype == 'potential'): # normalize around 0
      fmax = np.max(pool.map(maxabs_of_fieldblob, fnames))
      fmin = -fmax
    else:
      print_usage()
      sys.exit(1)

    print "fmin: {fmin}".format(fmin = fmin)
    print "fmax: {fmax}".format(fmax = fmax)

    pool.map(plot_fieldblob, ( (fn, fmin, fmax, ftype, cmname) for fn in fnames))
