#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

#from pylab import *
from scipy import optimize
from scipy import odr
import numpy as np
from numpy import *
import matplotlib.cm as cm
from matplotlib import pyplot as plt
from matplotlib.pyplot import *
import sys
import os as os
import json

print "Usage: %s [freqs] [noscale] [nofit] [showsep] [raitza]" % sys.argv[0]

freqsonly = ('freqs'   in sys.argv)
noscale   = ('noscale' in sys.argv)
nofit     = ('nofit'   in sys.argv)
showsep   = ('showsep' in sys.argv)
raitza    = ('raitza'  in sys.argv)

if (raitza):
  if (not os.path.exists('./raitza_ReKt.dat')):
    print "./raitza_ReKt.dat is missing - cannot plot his data"
    raitza = False

scalefac = 2.

f = open('./simparams.dat', 'r')
simparams = json.load(f)
f.close()

tau = 1.97867 # 95% confidence interval for 128 samples
              # Mathematica:
              # Needs["HypothesisTesting`"]; StudentTCI[0, 1, 128]

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

def enlargeticks():
	setp(gca().get_ymajorticklabels(), fontsize='large') 
	setp(gca().get_xmajorticklabels(), fontsize='large') 
	for l in gca().get_xticklines() + gca().get_yticklines(): 
	    l.set_markersize(6) 
	    l.set_markeredgewidth(1.2) 
	for l in gca().yaxis.get_minorticklines()+gca().xaxis.get_minorticklines():
	    l.set_markersize(3) 
	    l.set_markeredgewidth(1.2) 

#def myfunc(p, w):
#	return p[0] / (1 + (  (w**2-(p[2])**2)/(w*p[1])  )**2 )

def myfunc(p, w):
	return p[0] / (1 + ( (w-p[2]) / p[1] )**2 )**(3./2.)

# Parameters:
#   params[k,0] = K/nu
#   params[k,1] = nu
#   params[k,2] = omega0
def fitfunc(params, w):
	res = np.zeros(len(w))
        # make a rectangular array from linear param storage
	p   = np.reshape(params, (-1,3))

	# stupid way of implementing constraints
	if (any(p < 0.)):
		return res - 1.e100

	for k in range(len(p[:,0])):
		res = res + myfunc(p[k,:],w)

	return log(res)

def plotfitfuncs(params, xvals, scalefac):
        # make a rectangular array from linear param storage
	p   = np.reshape(params, (-1,3))

	for k in range(len(p[:,0])):
	        plot(xvals/scalefac, myfunc(p[k,:], xvals), "-",  color='0.3', linewidth=1.5, alpha=0.65)

def do_fit(filename, col, lims, p0):
  raw    = loadtxt("%s_jackknife.dat"        % filename)
  rawstd = loadtxt("%s_jackknife_stddev.dat" % filename)
  rawx  = raw[   :,  0]
  rawy  = raw[   :,col]
  raws  = rawstd[:,col]

  xlims = lims[0]/scalefac
  
  datamask  = (rawx < xlims[0]*scalefac) | (rawx > xlims[1]*scalefac)
  rawx     = ma.compressed(ma.array(rawx, mask=datamask))
  rawy     = ma.compressed(ma.array(rawy, mask=datamask))
  raws     = ma.compressed(ma.array(raws, mask=datamask))    

  # actually perform the fitting process
  d = odr.RealData(rawx, log(rawy), sy=log(raws/rawx))
  m = odr.Model(fitfunc)
  o = odr.ODR(d, m, p0, maxit=10000)
  o.set_job(fit_type=2)
  o.set_iprint(init=2, iter=2, final=2)
  out = o.run()

  p1       = out.beta
  fitstdev = out.sd_beta

  fitresults = np.reshape(p1, (-1, 3))
  fitstandarddevs = np.reshape(fitstdev, (-1, 3))

  # see http://docs.scipy.org/doc/scipy/reference/generated/scipy.odr.Output.html#scipy.odr.Output
  print "Initial values:"
  print np.reshape(p0, (-1,3))
  print 'Fitted parameters:'
  print fitresults
  print 'Scaled error bars:'
  print fitstandarddevs
  print 'Sum of squares error:'
  print out.sum_square
  #print 'Residual variance:'
  #print out.res_var
  print 'Relative error in function values computed within fcn.:'
  print out.rel_error
  print 'Stopreason:'
  print out.stopreason

  # calculate chi^2 ??
  #print np.sum(((dataY-fitfunc(p1, dataX))**2)/dataStDev)

  # pretty print important results
  #out.pprint()

  np.savetxt("fitresults.dat", fitresults)
  np.savetxt("fitstandarddevs.dat", fitstdev)

  return p1


def plotfile_momentum_electrons(filename, col, lims, ylabel, frequencies, params1, text1, params2, text2):
  raw    = loadtxt("%s_jackknife.dat"        % filename)
  rawstd = loadtxt("%s_jackknife_stddev.dat" % filename)
  numfreqs  = frequencies.size
  
  rawx  = raw[   :,  0]
  rawy  = raw[   :,col]
  raws  = rawstd[:,col]

  xlims = lims[0]/scalefac
  ylims = lims[1]
  
  if (not noscale):
    datamask  = (rawx < xlims[0]*scalefac) | (rawx > xlims[1]*scalefac)
    rawx     = ma.compressed(ma.array(rawx, mask=datamask))
    rawy     = ma.compressed(ma.array(rawy, mask=datamask))
    raws     = ma.compressed(ma.array(raws, mask=datamask))    
  
  rawyl = rawy - tau*raws
  rawyu = rawy + tau*raws

  fig = plt.figure(figsize=[11,8])

  if (raitza):
    raitzadata = np.loadtxt('./raitza_ReKt.dat')
    raitzawpl  = np.loadtxt('./raitza_wpl.dat')
    rxscale    = wpl/raitzawpl
    ryscale    = max(rawy)/max(raitzadata[:,1]) / 2.
    plot(raitzadata[:,0]*rxscale, raitzadata[:,1]*ryscale, "-", color='black',  label='Raitza', linewidth=2.0, alpha=1.0)
    
  # plot min-max value of confidence region
  fill_between(rawx/scalefac, rawyl.clip(min=ylims[0]), rawyu, color='blue', alpha=0.3)
  # plot raw data
  plot(rawx/scalefac, rawy, "o",  color='blue', label='raw data and 95% conf. int.', markersize=2)
    
  for ifreq in range(numfreqs):
    frequency = frequencies[ifreq]
    freqidx   = np.nonzero(rawx >= frequency)[0][0]
    realfreq  = rawx[freqidx]
    yval      = ylims[1]
    
    gca().annotate("",
            xy    = (realfreq, yval), xycoords='data',
            xytext= (realfreq, 0.2*yval), textcoords='data',
            arrowprops=dict(arrowstyle="->",connectionstyle="arc3",color='0.5'),
            )

  gca().set_yscale('log')
  gca().set_xscale('log')
  
  formatter = FuncFormatter(log_10_product)
  if not freqsonly:
    gca().xaxis.set_major_formatter(formatter)
    gca().set_yticklabels([])
  #gca().yaxis.set_major_formatter(formatter)

  gca().set_xticks(range(1,100))

  if (not noscale):
    gca().set_xlim(xlims)
    gca().set_ylim(ylims)
  
  gca().set_ylabel(ylabel, fontsize=25)
  gca().set_xlabel(r'$\omega\,[\mathrm{fs}^{-1}]$', fontsize=23)
  
  gca().plot([wpl,  wpl],  lims[1], '--', linewidth=2.0, color='#555555', alpha=0.75)
  gca().text(1.05*wpl, 1.075*lims[1][0], r'$\omega_\mathrm{pl}$', fontsize=23, horizontalalignment='left', verticalalignment='bottom', color='#555555')

  gca().plot([wMie,  wMie],  lims[1], '--', linewidth=2.0, color='#555555', alpha=0.75)
  gca().text(1.05*wMie, 1.075*lims[1][0], r'$\omega_\mathrm{Mie}$', fontsize=23, horizontalalignment='left', verticalalignment='bottom', color='#555555')

  if (showsep):
    plotfitfuncs(params2, rawx, scalefac)

  if (noscale):
    datamask  = (rawx < xlims[0]*scalefac) | (rawx > xlims[1]*scalefac)
    rawx     = ma.compressed(ma.array(rawx, mask=datamask))
    rawy     = ma.compressed(ma.array(rawy, mask=datamask))
    raws     = ma.compressed(ma.array(raws, mask=datamask))    

  # plot start parameters of fitting process
  if (nofit):
    plot(rawx/scalefac, exp(fitfunc(params1, rawx)), "-", color='green', label=text1, linewidth=2.5, alpha=0.65)
  else:
    plot(rawx/scalefac, exp(fitfunc(params2, rawx)), "-", color='red',  label=text2, linewidth=3.5, alpha=0.65)

  gca().legend(loc='lower left',scatterpoints=4,ncol=1,prop={'size':16},frameon=False)
  
  enlargeticks()

  if (not freqsonly):
    saveplot(filename)




#  initial guess for all of the peaks is read from file
if os.path.exists("./in_fitstartvalues.dat"):
  p0 = np.reshape(np.loadtxt("./in_fitstartvalues.dat", comments='#'),(-1))
else:
  p0 = [[1., 1., scalefac*wMie], [1., 1., scalefac*wpl]]
  np.savetxt("./in_fitstartvalues.dat", p0)
  

plotlims=np.loadtxt('./in_plotlims.dat')

if os.path.exists('./in_fitlims.dat'):
  fitlims=np.loadtxt('./in_fitlims.dat')
else:
  fitlims=plotlims
  np.savetxt('./in_fitlims.dat', fitlims)

if os.path.exists("./in_frequencies.dat"):
  frequencies=np.loadtxt("./in_frequencies.dat")
else:
  frequencies = []

if (nofit):
  p1 = p0
else:
  p1 = do_fit("momentum_electrons_acf", 4, fitlims, p0)

plotfile_momentum_electrons("momentum_electrons_acf", 4, [plotlims[0]*scalefac, fitlims[1]], "$\\log(\\mathrm{Re}\\{K(\\omega)\\})$", frequencies, p0, 'initial params', p1, "fit function")

plt.show()


