#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

from pylab import *
#from scipy import *
import numpy as numpy
from mpl_toolkits.mplot3d import axes3d

#from scipy import optimize
#from scipy import odr
import matplotlib.cm as cm
import sys

def saveplot(filenameout):
  plt.savefig("%s.pdf" % filenameout)
  plt.savefig("%s.png" % filenameout)


def plotfile_spherical_spectrum(filename, NR_NTheta_NPhi, lims):
  raw  = loadtxt("%s_jackknife.dat"        % filename)
  raws = loadtxt("%s_jackknife_stddev.dat" % filename) # TODO: plot stddev
  
  numcols = shape(raw)[1]
  
  xvals = raw[:,0]
  
  textx   = 0.80*lims[0][1]
  textidx = numpy.nonzero(xvals > textx)[0][0]
  textx   = xvals[textidx]
  
  fig = plt.figure(figsize=[15,10])
  
  fig.suptitle(filename)
  
  i = 0
  ntot = (NR_NTheta_NPhi[0]+1) * (NR_NTheta_NPhi[1]+1) * (NR_NTheta_NPhi[2]+1)
  
  for iR in range(0,NR_NTheta_NPhi[0]+1):
    for iTheta in range(0,NR_NTheta_NPhi[1]+1):
      for iPhi in range(0,NR_NTheta_NPhi[2]+1):
        i = i + 1

        plot(xvals, raw[:,i], "-",  color=cm.jet((1.*i)/ntot), linewidth=0.25, alpha=0.25)
   
        texty = raw[textidx,i]
    
        gca().text(textx, texty, "$iR= %d, i\Theta=%d, i\Phi=%d$" % (iR, iTheta, iPhi), color=cm.jet((1.*i)/ntot) )
    
  gca().set_xlim(lims[0])
  gca().set_ylim(lims[1])

  gca().set_yscale('log')
  gca().set_xscale('log')
  
  saveplot(filename)



def tbl(c, n, l, m):
  # symmetry relation, compare (7.23) in [Diss. Th. Raitza]
  if (m >= 0):
    return c[n-1,l,m].real
  else:
    return -c[n-1,l,abs(m)].conjugate().real
  
  

def plotfile_excitation_image_at_frequency_cone_theta(frequencies, filename, NR_NTheta_NPhi):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2.0
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  numfreqs  = len(frequencies)

  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  fig = plt.figure(figsize=(plotsize*numfreqs,NTheta*plotsize))
  fig.suptitle(filename)

  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  for ifreq in range(numfreqs):
    frequency = frequencies[ifreq]
    freqidx   = numpy.nonzero(xvals >= frequency)[0][0]
    realfreq  = xvals[freqidx]
    
    print "Processing frequency %d/%d: omega=%f fs^-1" % (ifreq+1, numfreqs, frequency)
    
    raw_spherical = np.empty((NR+1, NTheta+1, NPhi+1))
    
    i = 0
  
    for iR in range(0,NR_NTheta_NPhi[0]+1):
      for iTheta in range(0,NR_NTheta_NPhi[1]+1):
        for iPhi in range(0,NR_NTheta_NPhi[2]+1):
	  i = i + 1
          raw_spherical[iR, iTheta, iPhi] = raw[freqidx,i]
	  
    rvals     = linspace(0.,     1., num=NR    +1, endpoint=True)
    thetavals = linspace(0.,  np.pi, num=NTheta+1, endpoint=True)
    phivals   = linspace(0.,2*np.pi, num=NPhi  +1, endpoint=True)
    
    phiv, rv = meshgrid(phivals, rvals)
    phiv     = phiv.flatten()
    rv       =   rv.flatten()
    
    for iTheta in range(0,NTheta+1):
      ax = plt.subplot(NTheta+1, numfreqs, ifreq + 1 + numfreqs*iTheta, polar=True)
      ax.scatter(phiv, rv*sin(thetavals[iTheta]), c=raw_spherical[:,iTheta,:].flatten(), linewidth=0.)
        
      ax.set_xlabel("")#$\Phi$")
      ax.set_ylabel("")#$r$")
      ax.set_yticklabels([])
      ax.set_xticklabels([])
      ax.set_ylim([0.,1.])
      ax.set_title("$\Theta=%4.2f\pi$" % (thetavals[iTheta]/np.pi), fontsize=8)
      ax.set_alpha(0.25)
      if (iTheta==0):
        ax.set_title("$\omega=%6.4f\,\mathrm{fs}^{-1}$" % realfreq)

  saveplot("%s_excitation_at_frequency_coneTheta" % filename)


def plotfile_excitation_image_at_frequency_slice_phi(frequencies, filename, NR_NTheta_NPhi):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2.0
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  numfreqs  = len(frequencies)

  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  if ( (NPhi % 2) > 0 ):
    print "plotfile_excitation_image_at_frequency_slice_phi only supports even NPhi, but NPhi=%d" % NPhi
    return
  

  fig = plt.figure(figsize=(plotsize*numfreqs,NPhi/2*plotsize))
  fig.suptitle(filename)

  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  for ifreq in range(numfreqs):
    frequency = frequencies[ifreq]
    freqidx   = numpy.nonzero(xvals >= frequency)[0][0]
    realfreq  = xvals[freqidx]
    
    print "Processing frequency %d/%d: omega=%f fs^-1" % (ifreq+1, numfreqs, frequency)
    
    raw_spherical = np.empty((NR+1, NTheta+1, NPhi+1))
    
    i = 0
  
    for iR in range(0,NR_NTheta_NPhi[0]+1):
      for iTheta in range(0,NR_NTheta_NPhi[1]+1):
        for iPhi in range(0,NR_NTheta_NPhi[2]+1):
	  i = i + 1
          raw_spherical[iR, iTheta, iPhi] = raw[freqidx,i]
	  
    rvals     = linspace(0.,     1., num=NR    +1, endpoint=True)
    thetavals = linspace(0.,  np.pi, num=NTheta+1, endpoint=True)
    phivals   = linspace(0.,2*np.pi, num=NPhi  +1, endpoint=True)
    
    thetav, rv = meshgrid(thetavals, rvals)
    thetav     = thetav.flatten()
    rv         =   rv.flatten()
    
    for iPhi in range(0,NPhi/2+1):
      ax = plt.subplot(NPhi/2+1, numfreqs, ifreq + 1 + numfreqs*iPhi, polar=True)
                 # v--- we mirror and rotate theta here to get a better impression of a sphere-slice
      ax.scatter(-thetav+np.pi/2., rv, c=raw_spherical[:,:,iPhi       ].flatten(), linewidth=0.)
      ax.scatter( thetav+np.pi/2., rv, c=raw_spherical[:,:,iPhi+NPhi/2].flatten(), linewidth=0.)
        
      ax.set_xlabel("")#$\Theta$")
      ax.set_ylabel("")#$r$")
      ax.set_yticklabels([])
      ax.set_xticklabels([])
      ax.set_ylim([0.,1.])
      ax.set_title("$\Phi=%4.2f\pi$" % (phivals[iPhi]/np.pi), fontsize=8)
      ax.set_alpha(0.25)
      if (iPhi==0):
        ax.set_title("$\omega=%6.4f\,\mathrm{fs}^{-1}$" % realfreq)

  saveplot("%s_excitation_at_frequency_slicePhi" % filename)



def plotfile_excitation_image_at_frequency_3d(frequencies, filename, NR_NTheta_NPhi):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2.0
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  numfreqs  = len(frequencies)

  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  fig = plt.figure(figsize=(plotsize*numfreqs,2*plotsize))
  fig.suptitle(filename)

  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  for ifreq in range(numfreqs):
    frequency = frequencies[ifreq]
    freqidx   = numpy.nonzero(xvals >= frequency)[0][0]
    realfreq  = xvals[freqidx]
    
    print "Processing frequency %d/%d: omega=%f fs^-1" % (ifreq+1, numfreqs, frequency)
    
    x    = np.empty((NR+1)*(NTheta+1)*(NPhi+1) )
    y    = np.empty((NR+1)*(NTheta+1)*(NPhi+1) )
    z    = np.empty((NR+1)*(NTheta+1)*(NPhi+1) )
    vals = np.empty((NR+1)*(NTheta+1)*(NPhi+1) )
    
    i = -1
  
    for iR in range(0,NR+1):
      myR = 1./NR*iR
      for iTheta in range(0,NTheta+1):
        myTheta = np.pi/NTheta*iTheta
        for iPhi in range(0,NPhi+1):
	  myPhi = 2.*np.pi/NPhi*iPhi
	  
	  i = i + 1
	  
	  x[i]    = myR * cos(myPhi) * sin(myTheta)
	  y[i]    = myR * sin(myPhi) * sin(myTheta)
	  z[i]    = myR * cos(myTheta)
	  vals[i] = raw[freqidx,i+1]
	  
    ax = plt.subplot(1, numfreqs, ifreq + 1, projection='3d')
    ax.set_title("$\omega=%6.4f\,\mathrm{fs}^{-1}$" % realfreq)
    ax.scatter(x, y, z, c=vals, linewidth=0.)
        
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_ylabel("Z")
    ax.set_xticklabels([])
    ax.set_yticklabels([])
    ax.set_zticklabels([])
    ax.set_xlim([-1.,1.])
    ax.set_ylim([-1.,1.])
    ax.set_zlim([-1.,1.])
    #ax.set_title("$\Theta=%4.2f\pi$" % (thetavals[iTheta]/np.pi), fontsize=8)
    ax.set_alpha(0.25)

  saveplot("%s_excitation_at_frequency_3d" % filename)


NR_NTheta_NPhi = [16, 6, 6]	

# high-density case
#frequencies=[4.65, 6.48, 6.85, 7.18, 7.48, 7.81, 8.18]
# low-density case
frequencies=[4.33, 5.55, 6.37, 6.60, 6.68, 6.85, 7.13, 7.40]

plotfile_spherical_spectrum("field_spherical_spectrum", NR_NTheta_NPhi, [[3.0,15.0], [0.5, 300.0]])
plotfile_excitation_image_at_frequency_cone_theta(frequencies, "field_spherical_spectrum", NR_NTheta_NPhi)
plotfile_excitation_image_at_frequency_slice_phi(frequencies, "field_spherical_spectrum", NR_NTheta_NPhi)
plotfile_excitation_image_at_frequency_3d(frequencies, "field_spherical_spectrum", NR_NTheta_NPhi)



plt.show()
