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


def plotfile_spherical_fourier_coeffs(filename, NR_NTheta_NPhi, lims):
  raw  = loadtxt("%s_jackknife.dat"        % filename)
  raws = loadtxt("%s_jackknife_stddev.dat" % filename) # TODO: plot stddev
  
  numcols = shape(raw)[1]
  
  xvals = raw[:,0]
  
  textx   = 0.80*lims[0][1]
  textidx = numpy.nonzero(xvals > textx)[0][0]
  textx   = xvals[textidx]
  
  fig = plt.figure(figsize=[15,10])
  
  fig.suptitle(filename)
  
  n =  0
  l = -1
  m = -1
  
  maxl = NR_NTheta_NPhi[1]/2
  maxn = NR_NTheta_NPhi[0]/2
  
  for i in range(1, numcols):
    if (m==l):
      m = 0
      l = l+1
    else:
      m = m+1
      
    if (l > maxl):
      l = 0
      n = n + 1

    plot(xvals, raw[:,i], "-",  color=cm.prism(i), linewidth=0.5)
    
    texty = raw[textidx,i]
    
    gca().text(textx, texty, "n= %d, l=%d, m=%d" % (n, l, m), color=cm.prism(i) )
    
  gca().set_xlim(lims[0])
  gca().set_ylim(lims[1])

  gca().set_yscale('symlog', linthreshy=lims[2])
  gca().set_xscale('log')
  

def contourplotstuff(nrows, ncols, irow, icol, xvals, yvals, values, xlabel, ylabel):
    ax = plt.subplot(nrows, ncols, icol + 1 + ncols*(irow-1))
    ax.contour( xvals, yvals, values, 15, linewidths=0.5,colors='k')
    ax.contourf(xvals, yvals, values, 15, cmap=plt.cm.jet)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_yticklabels([])
    ax.set_xticklabels([])

  

def plotfile_excitation_image_at_frequency(frequencies, filename, NR_NTheta_NPhi, maxR, use_raitza_definition):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 3
  
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
  
    coeffs    = numpy.empty( (maxn+1, maxl+1, maxl+1) )

    for i in range(1, numcols):
      if (m==l):
        m = 0
        l = l+1
      else:
        m = m+1

      if (l > maxl):
        l = 0
        n = n + 1

      coeffs[n-1, l, m] = raw[freqidx,i]
      
    nx, ny = (100, 100)
    x      = numpy.linspace(-maxR, maxR, nx)
    y      = numpy.linspace(-maxR, maxR, ny)
    xv, yv = meshgrid(x, y)
  
    values = numpy.zeros((nx,ny))
  
    for i in range(nx):
      for j in range(ny):
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1): # TODO: m=-l..l
	      values[i,j] = values[i,j] + coeffs[n-1,l,m] * bfunc.bindings.my_bfunc(n, l, m, [xv[i,j], yv[i,j], 0.], use_raitza_definition, maxR)
	    
    contourplotstuff(6, numfreqs, 1, ifreq, xv, yv,     values,  'X', 'Y')
    contourplotstuff(6, numfreqs, 2, ifreq, xv, yv, abs(values), 'X', 'Y')

  
    values = numpy.zeros((nx,ny))
  
    for i in range(nx):
      for j in range(ny):
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1): # TODO: m=-l..l
	      values[i,j] = values[i,j] + coeffs[n-1,l,m] * bfunc.bindings.my_bfunc(n, l, m, [xv[i,j], 0., yv[i,j]], use_raitza_definition, maxR)
	    
    contourplotstuff(6, numfreqs, 3, ifreq, xv, yv,     values,  'X', 'Z')
    contourplotstuff(6, numfreqs, 4, ifreq, xv, yv, abs(values), 'X', 'Z')

  
    values = numpy.zeros((nx,ny))
  
    for i in range(nx):
      for j in range(ny):
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1): # TODO: m=-l..l
	      values[i,j] = values[i,j] + coeffs[n-1,l,m] * bfunc.bindings.my_bfunc(n, l, m, [0., xv[i,j], yv[i,j]], use_raitza_definition, maxR)
	    
    contourplotstuff(6, numfreqs, 5, ifreq, xv, yv,     values,  'Y', 'Z')
    contourplotstuff(6, numfreqs, 6, ifreq, xv, yv, abs(values), 'Y', 'Z')

  
	
NR_NTheta_NPhi = [16, 6, 6]	


#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_raitza_phi", NR_NTheta_NPhi, [[3.0,15.0], [-1.5e9, 1.5e10], 5.0e5])
plotfile_excitation_image_at_frequency([4.65, 6.48, 6.85, 7.18, 7.48, 7.81, 8.18], "field_spherical_fourier_coeffs_raitza_phi", NR_NTheta_NPhi, 1.0, True)

#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_wang_phi", NR_NTheta_NPhi,   [[3.0,15.0], [-8.0e-5, 5.0e-5 ], 3.0e-5])
#plotfile_excitation_image_at_frequency([4.65, 6.48, 6.85, 7.18, 7.48, 7.81, 8.18], "field_spherical_fourier_coeffs_wang_phi", NR_NTheta_NPhi, 1.0, False)

#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_raitza_ex", [[3.0,15.0], [-1.5e9, 1.5e10], 5.0e5])
#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_wang_ex",   [[3.0,15.0], [-1.5e5, 2.0e5 ], 2.5e2])

#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_raitza_ey", [[3.0,15.0], [-1.5e9, 1.5e10], 5.0e5])
#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_wang_ey",   [[3.0,15.0], [-1.5e5, 2.0e5 ], 2.5e2])

#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_raitza_ez", [[3.0,15.0], [-1.5e9, 1.5e10], 5.0e5])
#plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_wang_ez",   [[3.0,15.0], [-1.5e5, 2.0e5 ], 2.5e2])


plt.show()
