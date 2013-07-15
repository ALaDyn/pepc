#!/usr/bin/env python

import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if __name__ == '__main__':
  if (len(sys.argv) > 1):
    fname = sys.argv[1]
  else:
    fname = 'energy.csv'
  data = np.loadtxt(fname, skiprows = 1, delimiter = ',')
  t  = data[:,0]
  T  = data[:,1]
  Tc = data[:,2]
  V  = data[:,3]
  H  = data[:,4]

  fig = plt.figure()

  ax1 = fig.add_subplot(2, 1, 1)
  ax1.plot(t, H, 'k-', label = 'H')
  ax1.plot(t, T, 'k--', label = 'T')
  ax1.plot(t, np.max(H) * np.ones_like(t), 'k:', label = 'conserved')
  ax1.legend(loc = 'best')
  ax1.set_ylabel('E (a.u.)')

  ax2 = fig.add_subplot(2, 1, 2)
  ax2.plot(t, V, 'k-', label = 'V')
  ax2.plot(t, Tc, 'k--', label = 'Tc')
  ax2.legend(loc = 'best')
  ax2.set_ylabel('E (a.u.)')

  ax2.set_xlabel('omega_{p,e} t')

  fig.savefig('energy.png')
  fig.savefig('energy.pdf')
