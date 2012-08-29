#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

from pylab import *
#from scipy import *
import numpy as numpy

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
  
  

def plotfile_excitation_image_at_frequency(frequencies, filename, NR_NTheta_NPhi, maxR, use_raitza_definition):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 3
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  
  maxl      = NR_NTheta_NPhi[1]/2
  maxn      = NR_NTheta_NPhi[0]/2
  numfreqs  = len(frequencies)

  fig = plt.figure(figsize=(plotsize*numfreqs,2*plotsize))
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
  
    nx, nz = (100, 100)
    x      = numpy.linspace(-maxR, maxR, nx)
    z      = numpy.linspace(-maxR, maxR, nz)
    xv, zv = meshgrid(x, z)
  
    values = numpy.zeros((nx,nz))
  
    for i in range(nx):
      for j in range(nz):
        # treat xv[i,j], yv[i,j]
        for n in range(1,maxn+1):
          for l in range(0,maxl+1):
            for m in range(0,l+1):
	      values[i,j] = values[i,j] + coeffs[n-1,l,m] * bfunc.bindings.my_bfunc(n, l, m, [xv[i,j], 0., zv[i,j]], use_raitza_definition, maxR)
	    
    plt.subplot(2, numfreqs, ifreq + 1)
    plt.contour( xv, zv, values, 15, linewidths=0.5,colors='k')
    plt.contourf(xv, zv, values, 15, cmap=plt.cm.jet)
    plt.title("$\omega = %f fs^{-1}$" % realfreq)
  
    plt.subplot(2, numfreqs, ifreq + 1 + numfreqs)
    plt.contour( xv, zv, abs(values), 15, linewidths=0.5,colors='k')
    plt.contourf(xv, zv, abs(values), 15, cmap=plt.cm.jet)
  
  
  

  
	
NR_NTheta_NPhi = [16, 6, 6]	


plotfile_spherical_fourier_coeffs("field_spherical_fourier_coeffs_raitza_phi", NR_NTheta_NPhi, [[3.0,15.0], [-1.5e9, 1.5e10], 5.0e5])
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
