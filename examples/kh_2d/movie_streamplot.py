#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import fieldblob as fb

def plot_fieldblob2D(arg):
  fnames, fmin, fmax, cmname = arg
  print "plotting: %s vs. %s" % (fnames[0], fnames[1])
  dpi = 80.
  
  nx   = fb.n_of_fieldblob(fnames[0])
  tx   = fb.t_of_fieldblob(fnames[0])
  u   = fb.field_of_fieldblob(fnames[0])

  ny   = fb.n_of_fieldblob(fnames[1])
  ty   = fb.t_of_fieldblob(fnames[1])
  v   = fb.field_of_fieldblob(fnames[1])
  
  if (any(nx != ny) or tx != ty):
    print_usage()
    sys.exit(1)

  fig = plt.figure(figsize = (nx[0] / dpi, nx[1] / dpi))
  cm  = plt.get_cmap(cmname)

  ax  = fig.add_axes([0,0,1,1], frameon = False)
  
  offset = fb.offset_of_fieldblob(fnames[0])
  extent = fb.extent_of_fieldblob(fnames[0])
  
  X = np.linspace(offset[0], extent[0] + offset[0], num=nx[0])
  Y = np.linspace(offset[1], extent[1] + offset[1], num=nx[1])
 
  # http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.streamplot 
  ax.streamplot(Y, X, 
                u, v, 
                color = u,
                cmap  = cm,
                norm  = Normalize(vmin = fmin, vmax = fmax)
                )

  ax.set_xticks([])
  ax.set_yticks([])

  ax.text(0.05, 0.95, 'w_{{p,e}} t = {0:6.3f}'.format(tx), transform = ax.transAxes),

  fig.savefig(
    os.path.join(os.path.dirname(fnames[0]),
      'streamplot_' + os.path.basename(fnames[0]) + '_' + os.path.basename(fnames[1]) + '.png'),
    dpi = dpi
  )
  
  plt.show()

def maxabs_of_fieldblob(fname):
  return np.max(np.abs(fb.field_of_fieldblob(fname)))

def print_usage():
  print "Usage: movie_streamplot.py colormap fieldblob1X [fieldblob2X ...] fieldblob1Y [fieldblob2Y ...]"
  print ""
  print "colormap : any color map known to matplotlib, try 'Blues'/'Reds' or 'RdBu_r'"

if __name__ == '__main__':
  import multiprocessing as mp
  if (len(sys.argv) < 4 or len(sys.argv) % 2 == 1):
    print_usage()
    sys.exit(1)
  else:
    cmname  = sys.argv[1]
    fnames  = np.reshape(np.array(sys.argv[2:]), (2,-1)).transpose()

    pool    = mp.Pool()

    fmax = np.max(pool.map(maxabs_of_fieldblob, fnames[:,0]))
    fmin = -fmax

    print "fmin: {fmin}".format(fmin = fmin)
    print "fmax: {fmax}".format(fmax = fmax)

    pool.map(plot_fieldblob2D, ( (fn, fmin, fmax, cmname) for fn in fnames))
