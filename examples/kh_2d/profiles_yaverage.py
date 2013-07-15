#!/usr/bin/env python

import os, sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import fieldblob as fb
import pubmpl

def print_usage():
  print "Usage: nice_yaverage.py nt"

if __name__ == '__main__':
  if (len(sys.argv) < 2):
    print_usage()
    sys.exit(1)
  else:
    groups = { "density": ["ne", "ni"], "velocity": ["vey", "viy"], "elfield": ["ex"] }
    labels = { "ne": "electrons", "ni": "ions", 
      "vey": "electrons", "viy": "ions" }
    ylabels = { "density": r"$n_\sigma$", "velocity": r"$v_{y,\sigma} / v_{\mathrm{th},e}$",
      'elfield': r'$E_x$' }
    linestyle = { 'ne': 'k-', 'ni': 'k:', 'vey': 'k-', 'viy': 'k:', 'ex': 'k-' }

    nt = int(sys.argv[1])

    pubmpl.set_params(relwidth = 1.0, aspect_ratio = 1.0)
    fig, axs = plt.subplots(len(groups), sharex = True)
    axs[-1].set_xlabel(r'$x / \lambda_D$')

    for igroup, group in enumerate(['density', 'elfield', 'velocity']):
      axs[igroup].set_ylabel(ylabels[group])

      for member in groups[group]:
        fname = 'fields/' + member + '_' + '{0:06d}'.format(nt) + '.bin'
        n   = fb.n_of_fieldblob(fname)
        t   = fb.t_of_fieldblob(fname)
        y   = fb.field_of_fieldblob(fname)
        if (member == 'ne'): y *= np.abs(fb.qe_of_fieldblob(fname))
        elif (member == 'ni'): y *= np.abs(fb.qi_of_fieldblob(fname))

        ymean = np.mean(y, axis = 1)
        offset = fb.offset_of_fieldblob(fname)
        extent = fb.extent_of_fieldblob(fname)
        dx = extent / n

        x = np.linspace(
          offset[0] + 0.5 * dx[0],
          offset[0] + extent[0] - 0.5 * dx[0],
          num = ymean.shape[0]
        )

        if (member in labels.keys()):
          axs[igroup].plot(x, ymean, linestyle[member], label = labels[member])
        else:
          axs[igroup].plot(x, ymean, linestyle[member])
        axs[igroup].set_xlim(offset[0], offset[0] + extent[0])

      axs[igroup].legend(loc = 'best')

    axs[0].set_title(r'$\omega_{p,e} \times t = ' + '{0:3.2f}'.format(fb.t_of_fieldblob(fname)) + r'$')
    #plt.tight_layout(pad = 0.4, rect = [0.0, 0.0, 0.5, 1.0])
    plt.tight_layout(pad = 0.4)

    fig.savefig('profiles_' + '{0:06d}'.format(nt) + '.png')
    fig.savefig('profiles_' + '{0:06d}'.format(nt) + '.pdf')
