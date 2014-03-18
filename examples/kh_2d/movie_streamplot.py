#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Qt4Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import re as re

import fieldblob as fb

# arrows are not plotted at every grid point but approximately at every VERTEXSTEPth point (scaled for the directions)
VERTEXSTEP=8
# length scaling factor for the arrows
LENGTHSCALE=10


def plot_fieldblob2D(arg):
  fnames, fmin, fmax, cmname = arg
  print "plotting: %s vs. %s" % (fnames[0], fnames[1])
  dpi = 80.
  
  nx   = fb.n_of_fieldblob(fnames[0])
  tx   = fb.t_of_fieldblob(fnames[0])
  u    = np.transpose(fb.field_of_fieldblob(fnames[0]))

  ny   = fb.n_of_fieldblob(fnames[1])
  ty   = fb.t_of_fieldblob(fnames[1])
  v    = np.transpose(fb.field_of_fieldblob(fnames[1]))
  
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
  
  VERTEXSTEPX=VERTEXSTEP
  VERTEXSTEPY=(VERTEXSTEPX*nx[1])/nx[0]
 
  # streamplot is not supported until matplotlib version 1.3.1
  # and thus not available on Juvis and Juropa
  # Furthermore, it consumes a lot of memory and CPU power...
  # http://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.streamplot 
  #ax.streamplot(X, Y, 
  #              u, v, 
  #              color = u,
  #              cmap  = cm,
  #              norm  = Normalize(vmin = fmin, vmax = fmax)
  #              )
  
  veclen = np.sqrt(u[::VERTEXSTEPY,::VERTEXSTEPX]**2 + v[::VERTEXSTEPY,::VERTEXSTEPX]**2)
  
  ax.quiver(X[::VERTEXSTEPX],
            Y[::VERTEXSTEPY],
            u[::VERTEXSTEPY,::VERTEXSTEPX],
            v[::VERTEXSTEPY,::VERTEXSTEPX],
            veclen, # color
            angles='xy',
            scale_units='xy',
            scale = 1./LENGTHSCALE,
            width = 0.005,
            pivot = 'mid',
            cmap  = cm,
            norm  = Normalize(vmin = fmin, vmax = fmax)
            )

  ax.set_xticks([])
  ax.set_yticks([])

  ax.text(0.05, 0.95, 'w_{{p,e}} t = {0:6.3f}'.format(tx), transform = ax.transAxes),
  
  m = re.search(r"([a-zA-Z]+)_(\d{6})\.(\w*)", fnames[0])
  field0    = m.group(1)
  timefield = m.group(2)
  suffix    = m.group(3)
  m = re.search(r"([a-zA-Z]+)_(\d{6})\.(\w*)", fnames[1])
  field1    = m.group(1)
  
  fig.savefig(
    os.path.join(os.path.dirname(fnames[0]),
      'streamplot_%s-%s_%s.%s.png' % (field0, field1, timefield, suffix)),
    dpi = dpi
  )
  
def minmax_of_fieldblob(fname):
  veclensq = fb.field_of_fieldblob(fname[0])**2 + fb.field_of_fieldblob(fname[1])**2
  return [np.max(veclensq), np.min(veclensq)]

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
    minmax  = np.array(pool.map(minmax_of_fieldblob, fnames))
    fmax = np.sqrt(np.max(minmax[:,0]))
    fmin = np.sqrt(np.min(minmax[:,1]))

    print "fmin: {fmin}".format(fmin = fmin)
    print "fmax: {fmax}".format(fmax = fmax)

    pool.map(plot_fieldblob2D, ( (fn, fmin, fmax, cmname) for fn in fnames))
