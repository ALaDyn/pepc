#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

from pylab import *
#from scipy import *
import numpy as numpy
from mpl_toolkits.mplot3d import axes3d

#from scipy import optimize
#from scipy import odr
import matplotlib.cm as cm
import sys
import os as os
from matplotlib.colors import LinearSegmentedColormap
import json

print "Usage: %s [freqs]" % sys.argv[0]

freqsonly  = False
for arg in sys.argv:
	if (arg == 'freqs'):
		freqsonly = True


f = open('./simparams.dat', 'r')
simparams = json.load(f)
f.close()

wpl  = simparams['wpl']
wMie = simparams['wMie']

# see http://www.scipy.org/Cookbook/Matplotlib/CustomLogLabels
def log_10_product(x, pos):
     """The two args are the value and tick position.
     Label ticks with the product of the exponentiation"""
     return '%1i' % (x)

def saveplot(filenameout):
  plt.savefig("%s.pdf" % filenameout)
  plt.savefig("%s.png" % filenameout)


def plotfile_spherical_spectrum(filename, NR_NTheta_NPhi, lims, ylabel, frequencies):
  raw  = loadtxt("%s_jackknife.dat"        % filename)
  raws = loadtxt("%s_jackknife_stddev.dat" % filename) # TODO: plot stddev
  
  numcols = shape(raw)[1]
  numfreqs  = frequencies.size
  
  xvals = raw[:,0]
  
  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  fig = plt.figure(figsize=[11,8])
  
  # Using contourf to provide my colorbar info, then clearing the figure
  Z = [[0,0],[0,0]]
  levels = np.linspace(0.,1.)
  CS3 = plt.contourf(Z, levels, cmap=cm.jet)
  plt.clf()
  
  
  #fig.suptitle(filename)
  
  i = 0
  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  for iR in range(0,NR+1):
  
    mydata = np.zeros(raw[:,i].shape)
    numdat = 0
    
    for iTheta in range(0,NTheta+1):
      for iPhi in range(0,NPhi+1):
        i = i + 1
	numdat = numdat + 1
	mydata = mydata + raw[:,i]
	

    plot(xvals, mydata/numdat, "-",  color=cm.jet((1.*i)/ntot), linewidth=1.5)
    
  ylims = lims[1]

  for ifreq in range(numfreqs):
    frequency = frequencies[ifreq]
    freqidx   = numpy.nonzero(xvals >= frequency)[0][0]
    realfreq  = xvals[freqidx]
    yval      = ylims[1]
    
    print "Processing frequency %d/%d: omega=%f fs^-1" % (ifreq+1, numfreqs, frequency)
    print realfreq, yval, ylims
    
    gca().annotate("",
            xy    = (realfreq, yval), xycoords='data',
            xytext= (realfreq, 0.75*yval), textcoords='data',
            arrowprops=dict(color="0.5", arrowstyle="->",
                            connectionstyle="arc3"),
            )

  gca().set_yscale('log')
  gca().set_xscale('log')
  
  formatter = FuncFormatter(log_10_product)
  if not freqsonly:
    gca().xaxis.set_major_formatter(formatter)
  #gca().yaxis.set_major_formatter(formatter)

    gca().set_yticklabels([])
    gca().set_xticks(range(1,100))

  gca().set_xlim(lims[0])
  gca().set_ylim(lims[1])
  
  gca().set_ylabel(ylabel, fontsize=25)
  gca().set_xlabel(r'$\omega\,[\mathrm{fs}^{-1}]$', fontsize=23)
  
  gca().plot([wpl,  wpl],  lims[1], '--', linewidth=2.0, color='#555555', alpha=0.75)
  gca().text(1.05*wpl, 1.075*lims[1][0], r'$\omega_\mathrm{pl}$', fontsize=23, horizontalalignment='left', verticalalignment='bottom', color='#555555')

  gca().plot([wMie,  wMie],  lims[1], '--', linewidth=2.0, color='#555555', alpha=0.75)
  gca().text(1.05*wMie, 1.075*lims[1][0], r'$\omega_\mathrm{Mie}$', fontsize=23, horizontalalignment='left', verticalalignment='bottom', color='#555555')
  
  cax = plt.axes([0.165, 0.225, 0.175, 0.025])
  cbar = gcf().colorbar(CS3, ticks=[0., 1.], orientation='horizontal', cax=cax)
  cbar.ax.set_xticklabels([r'$\left|\vec{r}\right|=0$', r'$\left|\vec{r}\right|=r_\mathrm{cluster}$'], fontsize=18)# horizontal colorbar

  if (not freqsonly):
    saveplot(filename)


def plotfile_excitation_image_at_frequency_slice_phi(frequencies, filename, NR_NTheta_NPhi):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2.0
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  numfreqs  = frequencies.size

  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  if ( (NPhi % 2) > 0 ):
    print "plotfile_excitation_image_at_frequency_slice_phi only supports even NPhi, but NPhi=%d" % NPhi
    return
  

  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  for ifreq in range(numfreqs):
    fig = plt.figure(figsize=(plotsize,plotsize))
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
    rv         =     rv.flatten()
    
    for iPhi in range(0,NPhi/2+1):
      raw_spherical[:,:,0       ] = raw_spherical[:,:,0       ] + raw_spherical[:,:,iPhi       ]
      raw_spherical[:,:,0+NPhi/2] = raw_spherical[:,:,0+NPhi/2] + raw_spherical[:,:,iPhi+NPhi/2]
    
    ax = plt.subplot(1, 1, 1, polar=True)
               # v--- we mirror and rotate theta here to get a better impression of a sphere-slice
    ax.scatter(-thetav+np.pi/2., rv, c=raw_spherical[:,:,0       ].flatten(), linewidth=0.,cmap=cm.BuPu)
    ax.scatter( thetav+np.pi/2., rv, c=raw_spherical[:,:,0+NPhi/2].flatten(), linewidth=0.,cmap=cm.BuPu)
        
    ax.set_xlabel("")#$\Theta$")
    ax.set_ylabel("")#$r$")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    ax.set_ylim([0.,1.])
    ax.set_alpha(0.25)
    fig.suptitle("$\omega=%4.2f\,\mathrm{fs}^{-1}$" % realfreq)

    saveplot("%s_excitation_at_frequency_slicePhi_ifreq%d" % (filename, ifreq))
    
def plotfile_smooth_excitation_image_at_frequency_slice_phi(frequencies, filename, NR_NTheta_NPhi):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2.0
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  numfreqs  = frequencies.size

  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  if ( (NPhi % 2) > 0 ):
    print "plotfile_excitation_image_at_frequency_slice_phi only supports even NPhi, but NPhi=%d" % NPhi
    return
  

  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  for ifreq in range(numfreqs):
    fig = plt.figure(figsize=(plotsize,plotsize))
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
    
    for iPhi in range(0,NPhi/2+1):
      raw_spherical[:,:,0       ] = raw_spherical[:,:,0       ] + raw_spherical[:,:,iPhi       ]
      raw_spherical[:,:,0+NPhi/2] = raw_spherical[:,:,0+NPhi/2] + raw_spherical[:,:,iPhi+NPhi/2]
    
    ax = plt.subplot(1, 1, 1, polar=True)
               # v--- we mirror and rotate theta here to get a better impression of a sphere-slice
    ax.pcolormesh(-thetav+np.pi/2, rv, raw_spherical[:,:,0       ].flatten(),shading='gouraud', cmap=cm.BuPu)
    ax.pcolormesh( thetav+np.pi/2, rv, raw_spherical[:,:,0+NPhi/2].flatten(),shading='gouraud', cmap=cm.BuPu)
        
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    #ax.set_xlim([-1.,1.])
    #ax.set_ylim([-1.,1.])
    fig.suptitle("$\omega=%4.2f\,\mathrm{fs}^{-1}$" % realfreq)

    saveplot("%s_excitation_at_frequency_slicePhi_ifreq%d_smooth" % (filename, ifreq))    


def plotfile_excitation_image_at_frequency_radial_dependence(frequencies, filename, NR_NTheta_NPhi, mywpl):
  raw  = loadtxt("%s_jackknife.dat" % filename)
  
  plotsize = 2.0
  
  cdict  = {'red':   ((0.0, 0.0, 0.0),
                      (0.8, 1.0, 1.0),
                      (1.0, 0.5, 0.5)),

            'green': ((0.0, 0.0, 0.0),
                      (0.8, 0.0, 0.0),
                      (1.0, 0.7, 0.7)),
 
            'blue':  ((0.0, 0.0, 0.0),
                      (0.9, 0.0, 0.0),
                      (1.0, 0.3, 0.3))
        }  
  blue_red = LinearSegmentedColormap('BlueRed', cdict)
  
  numcols   = shape(raw)[1]
  xvals     = raw[:,0]
  numfreqs  = frequencies.size

  NR     = NR_NTheta_NPhi[0]
  NTheta = NR_NTheta_NPhi[1]
  NPhi   = NR_NTheta_NPhi[2]

  if ( (NPhi % 2) > 0 ):
    print "plotfile_excitation_image_at_frequency_radial_dependence only supports even NPhi, but NPhi=%d" % NPhi
    return
  
  markers=['o', 'v', '^', 'p', 's', '*', 'x', 'D', '+']
  
  ntot = (NR+1) * (NTheta+1) * (NPhi+1)
  
  fig = plt.figure(figsize=(12,8))

  # Using contourf to provide my colorbar info, then clearing the figure
  Z = [[0,0],[0,0]]
  levels = np.linspace(0.,1.)
  CS3 = plt.contourf(Z, levels, cmap=blue_red)
  plt.clf()

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
	  
    rvals     = linspace(0.,     1., num=NR    +1, endpoint=True) #* simparams['rion']

    mydata = np.zeros(shape(rvals))    
    
    
    for iTheta in range(0,NR_NTheta_NPhi[1]+1):
      for iPhi in range(0,NPhi/2+1):
        mydata = mydata + raw_spherical[:,iTheta,iPhi] + raw_spherical[:,iTheta,iPhi+NPhi/2]
    
    freqcolor = (realfreq/mywpl - sqrt(1./3.)) / (1-sqrt(1./3.))
    print freqcolor
    freqcolor=min(max(freqcolor, 0.), 1.)
    #gca().plot(rvals, mydata, label="$\omega=%4.2f\,\mathrm{fs}^{-1}=%4.2f\,\omega_\mathrm{pl}^{(\mathrm{cl})}$" % (realfreq, realfreq/mywpl),linewidth=2.0,marker=markers[ifreq],markersize=7,color=cm.jet(freqcolor))
    gca().plot(rvals, mydata, label="$\omega=%4.2f\,\mathrm{fs}^{-1}$" % realfreq,linewidth=2.0,marker=markers[ifreq],markersize=7,color=blue_red(freqcolor))
      
  gca().set_xlabel("$r / r_\mathrm{cluster}$", fontsize=22)
  gca().set_ylabel("$\\tilde{\Phi}(r,\omega)\,[\mathrm{a.u.}]$", fontsize=22)
  gca().set_yticklabels([])
  gca().set_xlim([0, 1.0])
  #gca().set_yscale('log')
  # ax.set_xticklabels([])
  # ax.set_ylim([0.,1.])
  # ax.set_alpha(0.25)
  gca().legend(ncol=1,prop={'size':16})
  cax = plt.axes([0.375, 0.85, 0.3, 0.025])
  cbar = gcf().colorbar(CS3, ticks=[0., 1.], orientation='horizontal', cax=cax)
  cbar.ax.set_xticklabels([r'$\omega_\mathrm{Mie}^{(\mathrm{cl})}$', r'$\omega_\mathrm{pl}^{(\mathrm{cl})}$'], fontsize=20)# horizontal colorbar

  saveplot("%s_excitation_at_frequency_radial_dependence" % filename)




NR_NTheta_NPhi = [16, 6, 6]

plotlims=np.loadtxt('./in_plotlims.dat')

if os.path.exists("./in_frequencies.dat"):
  frequencies=np.loadtxt("./in_frequencies.dat")
else:
  frequencies = []

plotfile_spherical_spectrum("field_spherical_spectrum_phi", NR_NTheta_NPhi, plotlims, r'$\tilde{\Phi}(\omega,\vec{r})\,[\mathrm{a.u.}]$', frequencies)

if os.path.exists("./in_frequencies.dat") and (not freqsonly):
  frequencies=np.loadtxt("./in_frequencies.dat")
  plotfile_excitation_image_at_frequency_slice_phi(frequencies, "field_spherical_spectrum_phi", NR_NTheta_NPhi)
  plotfile_smooth_excitation_image_at_frequency_slice_phi(frequencies, "field_spherical_spectrum_phi", NR_NTheta_NPhi)
  plotfile_excitation_image_at_frequency_radial_dependence(frequencies, "field_spherical_spectrum_phi", NR_NTheta_NPhi, frequencies[-1])

if (freqsonly):
  plt.show()
