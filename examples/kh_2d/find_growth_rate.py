#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt

import fieldblob as fb
from fit import lsfit

expmodel = lambda x, p: p[0] * np.exp(x * p[1]) + p[2]**2

def spectrum_of_fieldblob(fname):
  y  = fb.field_of_fieldblob(fname)
  yc = y[256 - 25: 256 + 25]
  
  return np.mean(np.abs(np.fft.rfft(yc, axis = 1))**2, axis = 0)

if __name__ == '__main__':
  if (len(sys.argv) < 4): sys.exit('syntax: find_growth_rate.py ntstart ntend step')
  ntstart = int(sys.argv[1])
  ntend   = int(sys.argv[2])
  step    = int(sys.argv[3])

  nts = range(ntstart, ntend + 1, step)
  files = [ 'fields/ey_{0:06d}.bin'.format(nt) for nt in nts ]
  spectra = [ spectrum_of_fieldblob(f) for f in files ]
  spectra = np.array(spectra)
  t = np.fromiter(( fb.t_of_fieldblob(f) for f in files ),
    dtype = np.float64, count = len(nts))

  ly = fb.extenty_of_fieldblob(files[0])
  ny = fb.ny_of_fieldblob(files[0])
  k  = 2 * np.pi * np.fft.fftfreq(ny, ly / ny)

  i0 = np.unravel_index(np.argmax(spectra), spectra.shape)[1]
  k0 = k[i0]
  dk = k[i0 + 1] - k[i0]
  y = spectra[:, i0]
  p = lsfit(t, y, expmodel, [ 1.0, 0.01, 0.0 ], maxfev = 10000)[0]
  gamma0 = p[1]
  print 'i0:    ', i0
  print 'k0:    ', k0, '+-', dk
  print 'gamma: ', p[1]
  plt.figure()
  plt.semilogy(t, y)
  plt.semilogy(t, expmodel(t, p))
  
  gamma = np.zeros(32)
  res = []
  for i in range(1,32):
    try:
      p = lsfit(t, spectra[:, i], expmodel, [ 1.0, 0.01, 0.0 ], maxfev = 10000)[0]
      if ((spectra[-1, i] - spectra[0, i]) > 100.0 * spectra[0, i]):
        print 'p[', i, ']: ', p
        gamma[i] = p[1]
        plt.figure()
        plt.title(str(i))
        plt.semilogy(t, spectra[:, i])
        plt.semilogy(t, expmodel(t, p))
        res.append((k[i], gamma[i]))
      else:
        print 'p[', i, ']: ', p, ' rejected!'
        gamma[i] = 0
    except RuntimeError:
      gamma[i] = 0

  plt.figure()
  plt.semilogy(k[1:32], gamma[1:32], 'k.')

  print '{', ', '.join(['{{{0}, {1}}}'.format(e[0], e[1]) for e in res]), '}'
  
  plt.show()

