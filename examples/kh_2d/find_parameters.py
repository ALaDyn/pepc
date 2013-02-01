#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt

import fieldblob as fb
from fit import lsfit

# fit models
softstep  = lambda x, p: p[0] * 0.5 * (np.tanh((x - p[1]) / p[2]) + 1)
softstepp = lambda x, p: p[0] * 0.5 * ( 1 - np.tanh((x - p[1]) / p[2])**2 ) / p[2]

peak  = lambda x, p: p[0] * (1 - np.tanh((x - p[1]) / p[2])**2)
peakp = lambda x, p: -2 * peak(x, p) * np.tanh((x - p[1]) / p[2]) / p[2]

th  = lambda x, p: np.tanh((x - p[0]) / p[1]) - 1
thp = lambda x, p: (1 - np.tanh((x - p[0]) / p[1])**2) / p[1]

def xaxis_of_fieldblob(fname):
  nx = fb.nx_of_fieldblob(fname)
  ox = fb.offsetx_of_fieldblob(fname)
  lx = fb.extentx_of_fieldblob(fname)
  return np.linspace(ox, ox + lx, nx)

def find_parameters_density(fname):
  q   = fb.qe_of_fieldblob(fname) if ('ne' in fname) else fb.qi_of_fieldblob(fname)
  rho = q * np.mean(fb.field_of_fieldblob(fname), axis = 1)
  x   = xaxis_of_fieldblob(fname)

  p0 = [ np.sign(q), 0.0, 5.0 ] 
  p = lsfit(x, rho, softstep, p0)[0]
  plt.figure()
  plt.plot(x, rho, 'kx')
  plt.plot(x, softstep(x, p), 'k-')
  plt.xlabel('x (electron Debye length)')
  plt.ylabel('n (electron Debye length)^(-3)')
  plt.title('density profile fit: ' + fname)
  return p

def find_parameters_velocity(fname):
  vy = np.mean(fb.field_of_fieldblob(fname), axis = 1)
  x  = xaxis_of_fieldblob(fname)

  i0 = np.argmax(np.abs(vy))
  p0 = [ vy[i0], x[i0], 1.0 ]
  p = lsfit(x, vy, peak, p0)[0]
  plt.figure()
  plt.plot(x, vy, 'kx')
  plt.plot(x, peak(x, p), 'k-')
  plt.xlabel('x (eletron Debye length)')
  plt.ylabel('vy (electron thermal velocity)')
  plt.title('velocity profile fit: ' + fname)
  return p

def find_parameters_velocity_th(fname):
  vy = np.mean(fb.field_of_fieldblob(fname), axis = 1)
  x = xaxis_of_fieldblob(fname)

  i0 = np.argmax(np.abs(vy))
  p0 = [ x[i0] + 10.0, 1.0 ]
  p = np.zeros(3)
  p[0] = -0.5 * vy[i0]
  p[1:] = lsfit(x[i0:], -2 * vy[i0:] / vy[i0], th, p0)[0]
  plt.figure()
  plt.plot(x, vy, 'kx')
  plt.plot(x, p[0] * th(x, p[1:]), 'k-')
  plt.xlabel('x (eletron Debye length)')
  plt.ylabel('vy (electron thermal velocity)')
  plt.title('velocity profile fit: ' + fname)
  return p

if __name__ == '__main__':
  if (len(sys.argv) < 2): sys.exit("syntax: find_parameters.py nt")
  nt = int(sys.argv[1])
  ntstr = '{0:06d}'.format(nt)

  file_pairs  = [ ('fields/n' + sp + '_' + ntstr + '.bin',
    'fields/v' + sp + 'y_' + ntstr + '.bin' ) for sp in ('e', 'i') ]
  
  for fp in file_pairs:
    print 'Parameters for: ' + fp[0]
    prho = find_parameters_density(fp[0])
    rho0, x0rho, Delta = prho
    print 'rho0:    ', rho0
    print 'x0:      ', x0rho
    print 'Delta:   ', Delta
    print "rho'0:   ", softstepp(x0rho, prho)

    print 'Parameters for: ' + fp[1]
    pv = find_parameters_velocity(fp[1])
    v0, x0v, delta = pv
    print 'v0: ', peak(x0rho, pv)
    print 'delta:', delta
    print "v'0: ", peakp(x0rho, pv)

    print 'Teilhaber profile parameters:'
    pvi = find_parameters_velocity_th(fp[1])
    v0, x0v, delta = pvi
    print 'v0: ', v0
    print 'delta:', delta
    print "v'0: ", pvi[0] * thp(x0v, pvi[1:])

  plt.show()
