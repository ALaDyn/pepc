#!/usr/bin/env python

import sys
import glob
import multiprocessing as mp

from math import ceil

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy.stats import t as student_tau

import fieldblob as fb
from fit import lsfit

expmodel = lambda x, p: p[0] * np.exp(x * p[1]) + p[2]**2
linmodel = lambda x, p: np.log(p[0]) + x * p[1]

def ky_of_fieldblob(fn):
  ly = fb.extenty_of_fieldblob(fn)
  ny = fb.ny_of_fieldblob(fn)
  return 2 * np.pi * np.fft.fftfreq(ny, ly / ny)

def spectrum_of_fieldblob(fname):
  y  = fb.field_of_fieldblob(fname)
  yc = y[127:128]

  #return np.mean(np.abs(np.fft.rfft(yc, axis = 1))**2, axis = 0)
  return np.mean(np.abs(np.fft.rfft(yc, axis = 1)), axis = 0)

if __name__ == '__main__':
  if (len(sys.argv) < 3): sys.exit('syntax: find_growth_rate.py { potential | ey | ... } nmodes [tstart tend]')
  field = sys.argv[1]
  nmodes = int(sys.argv[2])
  files = glob.glob('fields/' + field + '_??????.bin')
  files.sort()

  pool = mp.Pool()

  t = np.array(pool.map(fb.t_of_fieldblob, files))
  (tstart, tend) = (np.float(sys.argv[3]), np.float(sys.argv[4])) if (len(sys.argv) > 3) else (t[0], t[-1])

  # TODO: read from configuration
  B0 = 2.0
  M  = 16.0
  T  = 1.0
  Dx = np.sqrt(M * T) / B0
  V0 = Dx * B0 / M

  k = Dx * ky_of_fieldblob(files[0])
  spectra = np.array(pool.map(spectrum_of_fieldblob, files)).transpose()

  mask = np.where(np.logical_and(t >= tstart, t <= tend))
  tfit = t[mask]

  CONFIDENCE_LEVEL = 0.95
  #confidence_scale = student_tau.interval(CONFIDENCE_LEVEL, len(tfit), loc=0, scale=1.0)[1]
  confidence_scale = 1.95

  gamma = np.zeros(nmodes + 1)
  stddev = np.zeros_like(gamma)
  res = []

  mpl.rcParams['xtick.labelsize'] = 'small'
  mpl.rcParams['ytick.labelsize'] = 'small'
  mpl.rcParams['axes.labelsize'] = 'small'

  fig = plt.figure()
  for i in range(1, nmodes + 1):
    ax = fig.add_subplot(int(ceil(nmodes / 3.0)), 3, i)
    ax.set_title("mode number: " + str(i), fontsize = 'small')
    ax.semilogy(t, spectra[i, :], 'k-')
    try:
      (p, covariance) = lsfit(tfit, np.log(spectra[i][mask]), linmodel, [ 1.0, 0.01 ], maxfev = 10000)
      ax.semilogy(tfit, np.exp(linmodel(tfit, p)), 'r--')
      if True: # ((spectra[-1, i] - spectra[0, i]) > 1.0 * spectra[0, i]):
        print 'p[', i, ']: ', p
        gamma[i] = Dx / V0 * p[1]
        stddev[i] = Dx / V0 * np.sqrt(covariance[1,1])
        res.append((k[i], gamma[i]))
      else:
        print 'p[', i, ']: ', p, ' rejected!'
        gamma[i] = 0
    except RuntimeError:
      gamma[i] = 0

  #plt.tight_layout()
  plt.savefig('modes_%s_overview.png' % field)
  plt.savefig('modes_%s_overview.pdf' % field)

  plt.figure()
  x_analytic, y_analytic = np.loadtxt('growth_rate_analytic.csv', delimiter = ',', dtype = np.float64, unpack = True)
  plt.plot(x_analytic, y_analytic, 'k-')

  x_kai, y_kai, ub_kai = np.loadtxt('growth_rate_kai_ex1.csv', delimiter = ',', dtype = np.float64, unpack = True)
  plt.errorbar(x_kai, y_kai, yerr = ub_kai - y_kai, fmt = '.', color = 'gray')

  plt.errorbar(k[1:len(gamma)], gamma[1:], yerr = confidence_scale * stddev[1:], fmt = 'k.')

  print '{', ', '.join(['{{{0}, {1}}}'.format(e[0], e[1]) for e in res]), '}'

  plt.savefig('modes_%s_growthrate.png' % field)
  plt.savefig('modes_%s_growthrate.pdf' % field)


