#! /lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import locale
import os
import sys
import matplotlib.font_manager
import matplotlib.cm as cm
import re as re

print "Plots histogram data"
print "optimized for pepc-mw frontend"

print "Usage: %s [startstep]" % sys.argv[0]

locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

filename='hist/histogram_ve.%6.6d'

markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

datanames=[r'$v_x$', r'$v_y$', r'$v_z$', r'$|\vec{v}|$']

####################################

fig       = plt.figure(figsize=(12,8))
titletext = fig.suptitle(filename, fontsize=26)
#fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
prop      = matplotlib.font_manager.FontProperties(size=16)

pltlist     = list()
initialized = False
numcols     = 3 # set this value to 4 to include distribution function for magnitude of v, but: it is computed without drift correction in the code and hence cannot comply with the theoretical curve

if len(sys.argv[:]) > 1:
  currstep = int(sys.argv[1])
else:
  currstep = 1
  


maxx        = 1.
maxy        = 0.1
minx        = 0.
miny        = 0.

def readparams(filename):
    f=open(filename,'r')
    lines = f.readlines()
    
    res = {}
    res["state"]     = lines[0]
    res["minvals"]   = np.array(re.split('\s*', lines[1])[1:-1], dtype=float)
    res["maxvals"]   = np.array(re.split('\s*', lines[2])[1:-1], dtype=float)
    res["mean"]      = np.array(re.split('\s*', lines[3])[1:-1], dtype=float)
    res["numbins"]   = np.array(re.split('\s*', lines[4])[1:-1], dtype=int)
    res["sigma"]     = np.array(re.split('\s*', lines[5])[1:-1], dtype=float)
    res["nvals_tot"] = int(lines[6])
    res["mbins"]     = int(lines[7])
    res["ncols"]     = int(lines[8])
    res["binwidth"]  = np.array(re.split('\s*', lines[9])[1:-1], dtype=float)
        
    f.close()
    
    return res
    
    
def maxwell(minx, maxx, sigma, idx):

    if (idx < 3):
      xvals = np.linspace(minx,maxx)
      yvals = 1./(np.sqrt(2.*np.pi)*sigma)            * np.exp( -0.5*(xvals/sigma)**2 )
    else:
      xvals = np.linspace(0,maxx)
      yvals = np.sqrt(2./np.pi) / sigma**3 * xvals**2 * np.exp( -0.5*(xvals/sigma)**2 ) * xvals**2 * np.sqrt(2./(np.pi*sigma**3))

    return xvals, yvals
    

def initfig(*args):
    global numlinesold
    global pltlist
    global initialized
    global ax
    global minx
    global maxx
    
    if (initialized): # seltsamerweise wird die Funktion zweimal aufgerufen - das fangen wir hier ab
      return ax,
      
    initialized = True

    ax = fig.add_subplot(111)
    plt.xlabel("v")
    plt.ylabel("f(v)")
    plt.grid(True, which='both')
    #ax.set_yscale('log')
    plt.minorticks_on()

    data   = np.loadtxt(filename % currstep)
    params = readparams( "%s.params" % (filename % currstep))

    for colidx in range(0,numcols):
	pltlist.extend( ax.plot(data[:,2*colidx+1],data[:,2*colidx+2],label=datanames[colidx], linewidth=0.,color=cm.jet(1./numcols*colidx), markersize=10.0, markeredgewidth=1.0,marker=markers[colidx]))


    for colidx in range(0,numcols):
        xvals, yvals = maxwell(minx, maxx, params['sigma'][colidx], colidx)
	pltlist.extend( ax.plot(xvals,yvals, '--', label='Maxwell %s' % datanames[colidx], linewidth=2.5,color=cm.jet(1./numcols*colidx)))#, markersize=1.0, markeredgewidth=1.0,marker=markers[colidx])

    plt.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc=2, borderaxespad=0.)#, bbox_to_anchor=(1.05, 1))
    
    return ax,
    
    
def updatefig(*args):
    global currstep
    global pltlist
    global maxx
    global maxy
    global minx
    global miny
    global titletext
    
    if (os.path.exists(filename % (currstep+1))):

      currstep = currstep + 1
      data     = np.loadtxt(filename % currstep)
    
      for colidx in range(numcols):
        pltlist[colidx].set_data(data[:,2*colidx+1],data[:,2*colidx+2])
        
        localmaxx = data[:,2*colidx+1].max()
        localmaxy = data[:,2*colidx+2].max()
        localminx = data[:,2*colidx+1].min()
        localminy = data[:,2*colidx+2].min()
        
        maxx = max(localmaxx,maxx)
        maxy = max(localmaxy,maxy)
        minx = min(localminx,minx)
        miny = min(localminy,miny)

      params = readparams( "%s.params" % (filename % currstep))
      titletext.set_text("%6.6d  - %s" % (currstep, params['state']))

      for colidx in range(0,numcols):
          xvals, yvals = maxwell(minx, maxx, params['sigma'][colidx], colidx)
          pltlist[colidx+numcols].set_data(xvals, yvals)

      ax.set_xlim([minx,maxx])
      ax.set_ylim([miny,maxy])
    return tuple(pltlist)

####################################

ani = animation.FuncAnimation(fig, updatefig, interval=50, blit=False, init_func=initfig)
plt.show()



