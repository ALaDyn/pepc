#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

import numpy as np
import matplotlib.cm as cm
from matplotlib import pyplot as plt
import sys
import os as os
import json

scalefac = 2.

print "Usage: %s [showwidth] [raitza]" % sys.argv[0]

showwidth = ('showwidth' in sys.argv)
raitza    = ('raitza'  in sys.argv)

dirs     = [ '55.fields_resolution_lowdens',
             '300.fields_resolution_lowdens',
             '1k.fields_resolution_lowdens',
             '3k.fields_resolution_lowdens',
	     '30k.fields_resolution_lowdens',
	     '300k.fields_resolution_lowdens']

if (raitza):
  dataraitza = np.loadtxt('in_raitza_peaks.dat')


wplx = np.array([1., 1.e12])
wply = np.array([1., 1.])

def saveplot(filenameout):
  plt.savefig("%s.pdf" % filenameout)
  plt.savefig("%s.png" % filenameout)

def myarrow(x,y,wd):
        arrowlen        = 0.01
        arrowidthfactor = 0.025
	return [x,    (1-arrowidthfactor)*x, (1+arrowidthfactor)*x,    x,    x,    (1-arrowidthfactor)*x, (1+arrowidthfactor)*x,    x],  [        y+wd,         y+wd-arrowlen,         y+wd-arrowlen, y+wd, y-wd,            y-wd+arrowlen,         y-wd+arrowlen, y-wd]


fig = plt.figure(figsize=(8,5))
fig.subplots_adjust(bottom=0.15, right=0.875)
ax = fig.add_subplot(111)

plt.plot(wplx, wply,             '--', color='#ee8d18', linewidth=1.5)
plt.plot(wplx, wply/np.sqrt(3.), '--', color='#31b404', linewidth=1.5)
ax.text(1.1e7, 1.0,              '$\omega_{\\rm pl}$',           color='#ee8d18', size=28)
ax.text(1.1e7, 1.0/np.sqrt(3.0), '$\mathbf{\omega_{\\rm Mie}}$', color='#31b404', size=28)

#for i in range(1, 10):
#	wres = 1./np.sqrt(1.+(i+1.)/i)
#	plt.plot(wplx, wply*wres, '-', color='%f' % (1. - 1./(i+1)))

didlegend = False

nlist = []

for dirname in dirs:
	if (os.path.exists("./%s/fitresults.dat" % dirname)):
	        f = open('./%s/simparams.dat' % dirname, 'r')
        	simparams = json.load(f)
	        f.close()
		raw_peaks = np.loadtxt("./%s/fitresults.dat" % dirname)

		nparts   = simparams['Nion']
		wpl      = simparams['wpl']
		wMie     = simparams['wMie']
		
		nlist.append(nparts)

		for i in range((len(raw_peaks[:,0]))):
			cl = (0.,0.,1.,0.1)
			wd = raw_peaks[i,1]/wpl/scalefac
			w  = raw_peaks[i,2]/wpl/scalefac
		
			# plot each identified peak with its width
			if (showwidth):
				x, y = myarrow(nparts, w, wd/2)
				plt.plot(x, y, '-', color=cl, linewidth=1)
			
			if (didlegend):
			  plt.plot(nparts, w, 'o', markerfacecolor=cl)
			else:
			  plt.plot(nparts, w, 'o', markerfacecolor=cl, label='this work')
			  didlegend = True

if (raitza):
	plt.plot(dataraitza[:,0], dataraitza[:,1], 'd', markersize=10, markerfacecolor='none', markeredgecolor='black', markeredgewidth=1.0, label='Raitza')
	
	for nparts in dataraitza[:,0]:
	  nlist.append(nparts)

plt.gca().set_xscale('log')
plt.gca().set_yscale('log')
plt.gca().set_xlabel("$N_{\\mathrm{ion}}$", size=20)
plt.gca().set_ylabel("$\\frac{\\omega}{\\omega_{\\mathrm{pl}}}$", size=30, rotation='horizontal')
plt.legend(loc='lower right', frameon=True)
plt.xticks(nlist, map(int,nlist))
plt.minorticks_off()
plt.gca().set_yticks(np.array(range(1,11))/10.)
plt.gca().set_yticklabels(np.array(range(1,11))/10.)
plt.gca().set_xlim([10.,1.0e7])
plt.gca().set_ylim([  0.21,1.5])


saveplot("all_peaks.pdf")

plt.show()

