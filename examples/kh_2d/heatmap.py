#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

import fieldblob as fb
import pubmpl

def maxabs_of_fieldblob(fname):
  return np.max(np.abs(fb.field_of_fieldblob(fname)))

def print_usage():
  print "Usage: heatmap.py nt"

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print_usage()
    sys.exit(1)

  nt = int(sys.argv[1])

  pubmpl.set_params(relwidth = 1.0, aspect_ratio = 32.0 / 37.0)
  fig, ((axne, axni), (axpot, axey)) = plt.subplots(2, 2, sharex = True, sharey = True)

  for field in ('ne', 'ni', 'potential', 'ey'):
    fname = 'fields/' + field + '_' + '{0:06d}'.format(nt) + '.bin'
    y = fb.field_of_fieldblob(fname)
    if (field[0] == 'n'): # normalize starting at 0
      y = np.abs(y)
      fmax = np.max(y)
      fmin = 0.0
      if (field[1] == 'e'):
        cmname = 'Blues'
        ax = axne
      else:
        cmname = 'Reds'
        ax = axni
    else: # normalize around 0
      fmax = np.max(np.abs(y))
      fmin = -fmax
      cmname = 'RdBu_r'
      if (field == 'potential'):
        ax = axpot
      else:
        ax = axey

    cm = plt.get_cmap(cmname)

    offset = fb.offset_of_fieldblob(fname)
    extent = fb.extent_of_fieldblob(fname)
    ax.imshow(
      y.transpose(),
      interpolation = 'nearest',
      cmap = cm,
      norm = Normalize(vmin = fmin, vmax = fmax),
      extent = (offset[0], extent[0] + offset[0], offset[1], extent[1] + offset[1])
      #extent = (offset[1], extent[1] + offset[1], offset[0], extent[0] + offset[0])
    )

    if (field in ('potential', 'ey')):
      ax.set_xlabel(r'$x / \lambda_{D,e}$')
    if (field in ('ne', 'potential')):
      ax.set_ylabel(r'$y / \lambda_{D,e}$')

  plt.tight_layout(pad = 0.4, w_pad = 0.0)
  fig.savefig('heatmap_' + '{0:06d}'.format(nt) + '.png', dpi = 72.27)
  fig.savefig('heatmap_' + '{0:06d}'.format(nt) + '.pdf', dpi = 320)
