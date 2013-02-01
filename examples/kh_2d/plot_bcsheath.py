#!/usr/bin/env python

import glob
import multiprocessing as mp
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

import fieldblob as fb

def peak_vs_t_of_fieldblob(fname):
  print 'processing: ' + fname
  t = fb.t_of_fieldblob(fname)
  y = fb.field_of_fieldblob(fname)
  return(t, np.mean(y[-1, :]))

if __name__ == '__main__':
  nefnames = glob.glob('fields/ne_??????.bin')
  nifnames = glob.glob('fields/ni_??????.bin')

  pool = mp.Pool()
  
  tmp = pool.map(peak_vs_t_of_fieldblob, nefnames)
  tmp.sort()
  tne = [ x[0] for x in tmp ]
  nep = [ x[1] for x in tmp ]
  tme = pool.map(peak_vs_t_of_fieldblob, nifnames)
  tmp.sort()
  tni = [ x[0] for x in tmp ]
  nip = [ x[1] for x in tmp ]

  fig = plt.figure()
  ax = fig.add_subplot(111)
  ax.plot(tni, nip, 'k-', label = 'ni')
  ax.plot(tne, nep, 'r:', label = 'ne')

  ax.set_xlabel('w_{p,e} t')
  ax.set_ylabel('sheath density')

  ax.legend(loc = 'best')

  fig.savefig('bcsheath.png')
