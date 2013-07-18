#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import fieldblob as fb
import pubmpl

def print_usage():
  print 'Usage: nice_yaverage.py nt'

def calc_ylim(y):
  a = np.min(y)
  b = np.max(y)
  c = 0.1 * (b - a)
  return (a - c, b + c)

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print_usage()
    sys.exit(1)
  else:
    nt = int(sys.argv[1])

    pubmpl.set_params(relwidth = 1.0)
    #fig, (axrho, axe, axve, axvi) = plt.subplots(4, sharex = True)
    #axvi.set_xlabel(r'$x$')
    fig, (axrho, axe) = plt.subplots(2, sharex = True)
    axe.set_xlabel(r'$x$')

    ### plot 1) charge density ###
    fname = 'fields/ni_' + '{0:06d}'.format(nt) + '.bin'
    rhoi = fb.qi_of_fieldblob(fname) * fb.field_of_fieldblob(fname)
    fname = 'fields/ne_' + '{0:06d}'.format(nt) + '.bin'
    rhoe = fb.qe_of_fieldblob(fname) * fb.field_of_fieldblob(fname)
    rho = rhoi + rhoe

    rhomean = np.mean(rho, axis = 1)
    x = fb.xaxis_of_fieldblob(fname)

    axrho.set_title(r'$t = ' + '{0:3.2f}'.format(fb.t_of_fieldblob(fname)) + r'$')
    axrho.set_ylabel(r'$\rho$')
    axrho.set_ylim(calc_ylim(rhomean))
    axrho.get_yaxis().get_major_locator().set_params(nbins = 5)
    
    axrho.plot(x, rhomean, 'k-')

    ### plot 2) el. field ###
    fname = 'fields/ex_' + '{0:06d}'.format(nt) + '.bin'
    ex = fb.field_of_fieldblob(fname)
    exmean = np.mean(ex, axis = 1)

    axe.set_ylabel(r'$E_x$')
    axe.set_ylim(calc_ylim(exmean))
    axe.get_yaxis().get_major_locator().set_params(nbins = 5)

    axe.plot(x, exmean, 'k-')

    ### plot 3) vey ###
#    fname = 'fields/vey_' + '{0:06d}'.format(nt) + '.bin'
#    ve = fb.field_of_fieldblob(fname)
#    vemean = np.mean(ve, axis = 1)
#
#    axve.set_ylabel(r'$v_\mathrm{e}$')
#    axve.set_ylim(calc_ylim(vemean))
#    axve.get_yaxis().get_major_locator().set_params(nbins = 5)
#
#    axve.plot(x, vemean, 'k-')
#
    ### plot 4) vey ###
#    fname = 'fields/viy_' + '{0:06d}'.format(nt) + '.bin'
#    vi = fb.field_of_fieldblob(fname)
#    vimean = np.mean(vi, axis = 1)
#
#    axvi.set_xlim(np.min(x), np.max(x))
#    axvi.set_ylabel(r'$v_\mathrm{i}$')
#    axvi.set_ylim(calc_ylim(vimean))
#    axvi.get_yaxis().get_major_locator().set_params(nbins = 5)
#    
#    axvi.plot(x, vimean, 'k-')

    plt.tight_layout(pad = 0.0)

    fig.savefig('profiles_' + '{0:06d}'.format(nt) + '.png')
    fig.savefig('profiles_' + '{0:06d}'.format(nt) + '.pdf')
