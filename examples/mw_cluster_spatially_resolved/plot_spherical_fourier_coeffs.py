#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

from pylab import *
#from scipy import *
import numpy as numpy
from mpl_toolkits.mplot3d import axes3d

#from scipy import optimize
#from scipy import odr
import matplotlib.cm as cm
import sys

import bfunc # self-generated module using f2py

def saveplot(filenameout):
  plt.savefig("%s.pdf" % filenameout)
  plt.savefig("%s.png" % filenameout)


def plotfile_spherical_fourier_coeffs(filename, NR_NTheta_NPhi, lims):
  raw  = loadtxt("%s_jackknife.dat"        % filename)
  raws = loadtxt("%s_jackknife_stddev.dat" % filename) # TODO: plot stddev
  
  numcols = shape(raw)[1]
  
  xvals = raw[:,0]
  
  textx   = 0.75*lims[0][1]
  textidx = numpy.nonzero(xvals > textx)[0][0]
  textx   = xvals[textidx]
  
  fig = plt.figure(figsize=[15,10])
  
  fig.suptitle(filename)
  
  n =  1
  l = -1
  m = -1
  
  maxl = NR_NTheta_NPhi[1]/2
  maxn = NR_NTheta_NPhi[0]/2
  
  for i in range(1, numcols, 2):
    if (m==l):
      m = 0
      l = l+1
    else:
      m = m+1
      
    if (l > maxl):
      l = 0
      m = 0
      n = n + 1

    plot(xvals, raw[:,i], "-",  color=cm.prism(i), linewidth=0.5)
    
    texty = raw[textidx,i]
    
    gca().text(textx, texty, "Re, n= %d, l=%d, m=%d" % (n, l, m), color=cm.prism(i) )
    
    plot(xvals, raw[:,i+1], "-",  color=cm.prism(i), linewidth=0.5)
    
    texty = raw[textidx,i+1]
    
    gca().text(textx, texty, "Im, n= %d, l=%d, m=%d" % (n, l, m), color=cm.prism(i) )

  gca().set_xlim(lims[0])
  gca().set_ylim(lims[1])

  gca().set_yscale('symlog', linthreshy=lims[2])
  gca().set_xscale('log')

  saveplot(filename)



def contourplotstuff(nrows, ncols, irow, icol, xvals, yvals, values, xlabel, ylabel, freqval):
    ax = plt.subplot(nrows, ncols, icol + 1 + ncols*(irow-1))
    ax.imshow( values, cmap=cm.jet)
    #ax.contour( xvals, yvals, values, 15, linewidths=0.5,colors='k')
    #ax.contourf(xvals, yvals, values, 15, cmap=plt.cm.jet)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    if (irow==1):
      ax.set_title("$\omega=%6.4f\,\mathrm{fs}^{-1}$" % freqval)



def tbl(c, n, l, m):
  # symmetry relation, compare (7.23) in [Diss. Th. Raitza]
  if (m >= 0):
    return  c[n-1,l,m].real
  else:
    return -c[n-1,l,abs(m)].conjugate().real
  
  

def plotfile_excitation_image_at_frequency(frequencies, filename, NR_NTheta_NPhi, maxR):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2
  nx, ny   = (25, 25)
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  
  maxl      = NR_NTheta_NPhi[1]/2
  maxn      = NR_NTheta_NPhi[0]/2
  numfreqs  = len(frequencies)

  fig = plt.figure(figsize=(plotsize*numfreqs,6*plotsize))
  fig.suptitle(filename)

  for ifreq in range(numfreqs):
    frequency = frequencies[ifreq]
    freqidx   = numpy.nonzero(xvals >= frequency)[0][0]
    realfreq  = xvals[freqidx]
    
    print "Processing frequency %d/%d: omega=%f fs^-1" % (ifreq+1, numfreqs, frequency)
  
    n =  1
    l = -1
    m = -1
  
    coeffs    = numpy.empty( (maxn+1, maxl+1, maxl+1), dtype='complex' )

    for i in range(1, numcols, 2):
      if (m==l):
        m = 0
        l = l+1
      else:
        m = m+1

      if (l > maxl):
        l = 0
	m = 0
        n = n + 1

      coeffs[n-1, l, m] = complex(raw[freqidx,i], raw[freqidx,i+1])
      
    x      = numpy.linspace(-maxR, maxR, nx)
    y      = numpy.linspace(-maxR, maxR, ny)
    xv, yv = meshgrid(x, y)
  
    values = numpy.zeros((nx,ny), dtype = 'complex')
  
    for i in range(nx):
      for j in range(ny):
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1):
	      values[i,j] = values[i,j] + tbl(coeffs,n,l,m) * bfunc.bindings.my_bfunc(n, l, m, [xv[i,j], yv[i,j], 0.], maxR)
	    
    contourplotstuff(6, numfreqs, 1, ifreq, xv, yv,     values.real,  'X', 'Y', realfreq)
    contourplotstuff(6, numfreqs, 2, ifreq, xv, yv, abs(values.real), 'X', 'Y', realfreq)

  
    values = numpy.zeros((nx,ny), dtype = 'complex')
  
    for i in range(nx):
      for j in range(ny):
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1):
	      values[i,j] = values[i,j] + tbl(coeffs,n,l,m) * bfunc.bindings.my_bfunc(n, l, m, [xv[i,j], 0., yv[i,j]], maxR)
	    
    contourplotstuff(6, numfreqs, 3, ifreq, xv, yv,     values.real,  'X', 'Z', realfreq)
    contourplotstuff(6, numfreqs, 4, ifreq, xv, yv, abs(values.real), 'X', 'Z', realfreq)

  
    values = numpy.zeros((nx,ny), dtype = 'complex')
  
    for i in range(nx):
      for j in range(ny):
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1):
	      values[i,j] = values[i,j] + tbl(coeffs,n,l,m) * bfunc.bindings.my_bfunc(n, l, m, [0., xv[i,j], yv[i,j]], maxR)
	    
    contourplotstuff(6, numfreqs, 5, ifreq, xv, yv,     values.real,  'Y', 'Z', realfreq)
    contourplotstuff(6, numfreqs, 6, ifreq, xv, yv, abs(values.real), 'Y', 'Z', realfreq)

  saveplot("%s_excitation_at_frequency" % filename)
  
	
NR_NTheta_NPhi = [16, 6, 6]	

# high-density case
#frequencies=[4.65, 6.48, 6.85, 7.18, 7.48, 7.81, 8.18]
# low-density case
frequencies=[4.32, 6.40, 6.64, 6.84, 6.68, 7.13, 7.40]

plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_phi", NR_NTheta_NPhi,   [[3.0,15.0], [-1.0e-5, 1.0e-5 ], 1.0e-7])
plotfile_excitation_image_at_frequency(frequencies, "field_spherical_fourier_coeffs_phi", NR_NTheta_NPhi, 1.0)

plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_ex", NR_NTheta_NPhi,   [[3.0,15.0], [-5.0e-7, 5.0e-7 ], 1.0e-8])
plotfile_excitation_image_at_frequency(frequencies, "field_spherical_fourier_coeffs_ex", NR_NTheta_NPhi, 1.0)

#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_ey",   [[3.0,15.0], [-1.5e5, 2.0e5 ], 2.5e2])

#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_ez",   [[3.0,15.0], [-1.5e5, 2.0e5 ], 2.5e2])


plt.show()
