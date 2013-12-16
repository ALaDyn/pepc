#!/usr/bin/env python

import multiprocessing as mp
import sys
import glob
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import fieldblob as fb
from mpl_toolkits.mplot3d import Axes3D
from fit import lsfit

def ky_of_fieldblob(fn):
  ly = fb.extenty_of_fieldblob(fn)
  ny = fb.ny_of_fieldblob(fn)
  return 2 * np.pi * np.fft.fftfreq(ny, ly / ny)

def spectrum_of_fieldblob(fname):
  y  = fb.field_of_fieldblob(fname)
  nx = fb.nx_of_fieldblob(fname)
  yc = y[nx / 2 - nx / 20: nx / 2 + nx / 20]
  
  return np.mean(np.abs(np.fft.rfft(yc, axis = 1))**2, axis = 0)

expmodel = lambda x, p: p[0] * np.exp(x * p[1]) + p[2]**2

if __name__ == '__main__':
  if (len(sys.argv) < 3): sys.exit('syntax: plot_modes.py { potential | ey | ... } nmodes [tstart tend]')
  field = sys.argv[1]
  nmodes = int(sys.argv[2])
  fn = glob.glob('fields/' + field + '_??????.bin')
  fn.sort()

  pool = mp.Pool()

  t = np.array(pool.map(fb.t_of_fieldblob, fn))
  if (len(sys.argv) > 3):
    tstart = np.float(sys.argv[3])
    tend = np.float(sys.argv[4])
  else:
    tstart = t[0]
    tend = t[-1]
  k = ky_of_fieldblob(fn[0])
  xi = np.array(pool.map(spectrum_of_fieldblob, fn)).transpose()

  mpl.rcParams['xtick.labelsize'] = 'small'
  mpl.rcParams['ytick.labelsize'] = 'small'
  mpl.rcParams['axes.labelsize'] = 'small'

  fig = plt.figure()
  ax = fig.gca(projection = '3d')

  gamma = -1.0 * np.ones(nmodes + 1)
  mask = np.where(np.logical_and(t >= tstart, t <= tend))
  tfit = t[mask]
  naccepted = 0
  for ik in range(1, nmodes + 1):
    xifit = xi[ik][mask]
    ximin = xifit[0]
    ximax = xifit[-1]
    p = lsfit(tfit, xifit, expmodel, [ 1.0, 0.01, 0.0 ], maxfev = 10000)[0]
    if ((ximax - ximin) > 100.0 * ximin):
      print 'p[', ik, ']: ', p
      naccepted += 1
      gamma[ik] = p[1]
      ax.plot(k[ik] * np.ones_like(tfit), tfit, expmodel(tfit, p), 'k--', 
        linewidth = 1.5 * mpl.rcParams['lines.linewidth'])
    else:
      print 'p[', ik, ']: ', p, ' rejected!'

  cmap = mpl.cm.get_cmap('autumn')
  norm = mpl.colors.Normalize(vmin = 0, vmax = np.max(gamma))

  ax.plot([k[1], k[nmodes]], tstart * np.ones(2), np.zeros(2), 'r--')
  ax.plot([k[1], k[nmodes]], tend * np.ones(2), np.zeros(2), 'r--')

  for ik in range(1, nmodes + 1):
    if (gamma[ik] >= 0.0):
      ax.plot(k[ik] * np.ones_like(t), t, xi[ik], '-', color = cmap(norm(gamma[ik])))
    else:
      ax.plot(k[ik] * np.ones_like(t), t, xi[ik], 'k-')

  ax.set_xlabel('ky')
  ax.set_ylabel('t')
  ax.set_zlabel('|' + field + '(ky, t)|**2')

  cax, kw = mpl.colorbar.make_axes(ax)
  cb = mpl.colorbar.ColorbarBase(cax, norm = norm, cmap = cmap)
  cb.set_label('gamma')

  plt.savefig('modes_' + field + '.png')
  plt.savefig('modes_' + field + '.pdf')

  mpl.rcParams['text.usetex']=True
  mpl.rcParams['text.latex.unicode']=True

  plt.figure()
  cmap2 = mpl.cm.get_cmap('Set1')
  iaccepted = 0
  for ik in range(1, nmodes + 1):
    if (gamma[ik] >= 0.0):
      plt.semilogy(t, xi[ik], '-',
        color = cmap2(float(iaccepted) / float(naccepted)),
        label = "$m = {0}$".format(ik))
      iaccepted += 1

  plt.legend(loc = 'best', prop = { 'size': 'small' })
  plt.xlabel(r"$\omega_{p,e} \, t$")

  mpl.rcParams['text.usetex']=False
  mpl.rcParams['text.latex.unicode']=False

  plt.savefig('modes_log_' + field + '.png')
  plt.savefig('modes_log_' + field + '.pdf')

  plt.figure()
  plt.semilogy(k[np.where(gamma > 0.0)], gamma[np.where(gamma > 0.0)], 'k.')
  offset = mpl.transforms.offset_copy(plt.gca().transData, fig = plt.gcf(), y = 6.0, units = 'points')
  for n in np.arange(nmodes + 1)[np.where(gamma > 0.0)]:
    plt.text(k[n], gamma[n], str(n), transform = offset)
  plt.xlabel('ky')
  plt.ylabel('gamma')

  plt.savefig('growth_rates_' + field + '.png')
  plt.savefig('growth_rates_' + field + '.pdf')
