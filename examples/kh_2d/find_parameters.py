#!/usr/bin/env python

import sys

import numpy as np
import matplotlib.pyplot as plt

import fieldblob as fb
import pubmpl
from fit import lsfit

# fit models
softstep  = lambda x, p: p[0] * 0.5 * (np.tanh((x - p[1]) / p[2]) + 1)
softstepp = lambda x, p: p[0] * 0.5 * ( 1 - np.tanh((x - p[1]) / p[2])**2 ) / p[2]

peak  = lambda x, p: p[0] * (1 - np.tanh((x - p[1]) / p[2])**2)
peakp = lambda x, p: -2 * peak(x, p) * np.tanh((x - p[1]) / p[2]) / p[2]

th  = lambda x, p: np.tanh((x - p[0]) / p[1]) - 1
thp = lambda x, p: (1 - np.tanh((x - p[0]) / p[1])**2) / p[1]

def calc_ylim(y):
  a = np.min(y)
  b = np.max(y)
  c = 0.1 * (b - a)
  return (a - c, b + c)

def find_parameters_density(fname, ax = None):
  q   = fb.qe_of_fieldblob(fname) if ('ne' in fname) else fb.qi_of_fieldblob(fname)
  rho = q * np.mean(fb.field_of_fieldblob(fname), axis = 1)
  x   = fb.xaxis_of_fieldblob(fname)

  p0 = [ np.sign(q), 0.0, 5.0 ] 
  p = lsfit(x, rho, softstep, p0)[0]

  if (ax == None):
    plt.figure()
    plt.title('density profile fit: ' + fname + ', wp x t = ' + str(fb.t_of_fieldblob(fname)))
    ax = plt.gca()
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$\rho$')
  else:
    ax.get_yaxis().get_major_locator().set_params(nbins = 5)
    ax.set_ylim(calc_ylim(rho))

  ax.plot(x, rho, 'k.')
  ax.plot(x, softstep(x, p), 'r-')
  ax.set_xlim(np.min(x), np.max(x))

  return p

def find_parameters_velocity(fname, ax = None):
  vy = np.mean(fb.field_of_fieldblob(fname), axis = 1)
  x  = fb.xaxis_of_fieldblob(fname)

  i0 = np.argmax(np.abs(vy))
  p0 = [ vy[i0], x[i0], 1.0 ]
  p = lsfit(x, vy, peak, p0)[0]

  if (ax == None):
    plt.figure()
    plt.title('velocity profile fit: ' + fname + ', wp x t = ' + str(fb.t_of_fieldblob(fname)))
    ax = plt.gca()
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$v_\sigma$')
  else:
    ax.get_yaxis().get_major_locator().set_params(nbins = 5)
    ax.set_ylim(calc_ylim(vy))

  ax.plot(x, vy, 'k.')
  ax.plot(x, peak(x, p), 'r-')
  ax.set_xlim(np.min(x), np.max(x))

  return p

def find_parameters_velocity_th(fname, ax = None):
  vy = np.mean(fb.field_of_fieldblob(fname), axis = 1)
  x = fb.xaxis_of_fieldblob(fname)

  i0 = np.argmax(np.abs(vy))
  p0 = [ x[i0] + 10.0, 1.0 ]
  p = np.zeros(3)
  p[0] = -0.5 * vy[i0]
  p[1:] = lsfit(x[i0:], -2 * vy[i0:] / vy[i0], th, p0)[0]

  if (ax == None):
    plt.figure()
    plt.title('velocity profile fit: ' + fname + ', wp x t = ' + str(fb.t_of_fieldblob(fname)))
    ax = plt.gca()
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'$v_\sigma$')
  else:
    ax.get_yaxis().get_major_locator().set_params(nbins = 5)
    ax.set_ylim(calc_ylim(vy))

  ax.plot(x, vy, 'k.')
  ax.plot(x, p[0] * th(x, p[1:]), 'r-')
  ax.set_xlim(np.min(x), np.max(x))

  return p

if __name__ == '__main__':
  if (len(sys.argv) < 2): sys.exit("syntax: find_parameters.py nt")
  nt = int(sys.argv[1])
  ntstr = '{0:06d}'.format(nt)

  pubmpl.set_params(relwidth = 1.0)
  fig, axs = plt.subplots(2, 2, sharex = True)

  file_pairs  = [ ('fields/n' + sp + '_' + ntstr + '.bin',
    'fields/v' + sp + 'y_' + ntstr + '.bin' ) for sp in ('e', 'i') ]
  
  for ip, fp in enumerate(file_pairs):
    print 'Parameters for: ' + fp[0]
    prho = find_parameters_density(fp[0], axs[ip][0])
    rho0, x0rho, Delta = prho
    print 'rho0:    ', rho0
    print 'X0:      ', x0rho
    print 'Delta:   ', Delta
    print "rho'0:   ", softstepp(x0rho, prho)

#    print 'Parameters for: ' + fp[1]
#    pv = find_parameters_velocity(fp[1])
#    v0, x0v, delta = pv
#    print 'v0:    ', peak(x0rho, pv)
#    print 'x0:    ', x0v
#    print 'delta: ', delta
#    print "v'0:   ", peakp(x0rho, pv)

    print 'Teilhaber profile parameters:'
    pvi = find_parameters_velocity_th(fp[1], axs[ip][1])
    v0, x0v, delta = pvi
    print 'v0:    ', v0
    print 'x0:    ', x0v
    print 'delta: ', delta
    print "v'0:   ", pvi[0] * thp(x0v, pvi[1:])

  axs[0][0].set_ylabel(r'$\rho_\mathrm{e}$')

  axs[0][1].set_ylabel(r'$v_\mathrm{e}$')
  axs[0][1].yaxis.set_ticks_position("right")
  axs[0][1].yaxis.set_ticks_position("both")
  axs[0][1].yaxis.set_label_position("right")

  axs[1][0].set_xlabel(r'$x$')
  axs[1][0].set_ylabel(r'$\rho_\mathrm{i}$')

  axs[1][1].set_xlabel(r'$x$')
  axs[1][1].set_ylabel(r'$v_\mathrm{i}$')
  axs[1][1].yaxis.set_ticks_position("right")
  axs[1][1].yaxis.set_ticks_position("both")
  axs[1][1].yaxis.set_label_position("right")

  fig.text(0.5, 1.0,
    r'$t = ' + '{0:3.2f}'.format(fb.t_of_fieldblob(file_pairs[0][0])) + r'$',
    horizontalalignment = 'center', verticalalignment = 'top')
  plt.tight_layout(pad = 0.0, w_pad = 0.4, rect = (0, 0, 1, 0.92))
#  plt.show()
  fig.savefig('parameters_' + ntstr + '.pdf')
  fig.savefig('parameters_' + ntstr + '.png')
