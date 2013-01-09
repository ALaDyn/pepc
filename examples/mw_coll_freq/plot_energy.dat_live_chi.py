#! /lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import locale
import os
import sys
import matplotlib.font_manager
import matplotlib.cm as cm

print "Usage: %s [datafilename]" % sys.argv[0]

print "Plots data from energy.dat"
print "optimized for pepc-mw frontend"
print "Creating energy.dat.pdf in current directory"

locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

filename='energy.dat'

additionalfiles = sys.argv[1:]

print "Adding information from files ", additionalfiles

#markers=['+','*','.','1','2','3','4','<','>','D','H','^','_','d','h','o','p','s','v','x','|','None',' ','']
markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

numlinesold = 0

datanames=['$t [\mathrm{fs}]$',
           '$U_{\mathrm{pot}}^{\mathrm{(total)}}$',
           '$U_{\mathrm{pot}}^{\mathrm{(near field)}}$',
           '$U_{\mathrm{pot}}^{\mathrm{(far field)}}$',
           '$U_{\mathrm{kin}}^{\mathrm{(electrons)}}$',
           '$U_{\mathrm{kin}}^{\mathrm{(ions)}}$',
           '$U_{\mathrm{kin}}^{\mathrm{(total)}}$',
           '$U_{\mathrm{kin,w/o\ drift}}^{\mathrm{(electrons)}}$',
           '$U_{\mathrm{kin,w/o\ drift}}^{\mathrm{(ions)}}$',
           '$U^{\mathrm{(tot)}}$',
           '$T_{\mathrm{e}}^0$',
           '$T_{\mathrm{e}}^{\mathrm{uncor}}$',
           '$\chi_{\mathrm{e}}$',
           '$\Delta T_{\mathrm{e}}$',
           '$T_{\mathrm{i}}^0$',
           '$T_{\mathrm{i}}^{\mathrm{uncor}}$',
           '$\chi_{\mathrm{i}}$',
           '$\Delta T_{\mathrm{i}}$'
          ]


####################################

fig = plt.figure(figsize=(16,10))
fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
prop = matplotlib.font_manager.FontProperties(size=16)

pltlist = list()

initialized = False

def initfig(*args):
    global numlinesold
    global pltlist
    global initialized
    global ax
    
    if (initialized): # seltsamerweise wird die Funktion zweimal aufgerufen - das fangen wir hier ab
      return ax,
      
    initialized = True

    ax = fig.add_subplot(111)
    plt.xlabel(datanames[0])
    plt.ylabel("Energy [Ry]")
    plt.grid(True, which='both')
    #ax.set_yscale('log')
    plt.minorticks_on()

    data     = np.loadtxt(filename)
    numlines = data.shape[0]
    numcols  = data.shape[1]
    numlinesold = numlines

    for colidx in range(1,18):
	pltlist.extend(ax.plot(data[:,0],data[:,colidx],label=datanames[colidx], linewidth=1.0))#,color=cm.lines(colidx)))#, markersize=1.0, markeredgewidth=1.0,marker=markers[colidx])

    for f in additionalfiles:
	raw = loadtxt(f)
        ax.plot(raw[:,0],raw[:,1],label="\\verb+%s+" % f, linewidth=1.0)#,color=cm.hsv(24*colidx))#, markersize=1.0, markeredgewidth=1.0,marker=markers[colidx])

    plt.legend(numpoints=1, prop=prop, ncol=1, frameon=True, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
    return ax,
    
    
def updatefig(*args):
    global numlinesold
    global pltlist
    
    data     = np.loadtxt(filename)
    numlines = data.shape[0]
    numcols  = data.shape[1]
    
    if (numlines <> numlinesold):
      numlinesold = numlines
    
      for colidx in range(1,18):
        pltlist[colidx-1].set_data(data[:,0],data[:,colidx])

    ax.relim()
    ax.autoscale_view()
    return tuple(pltlist)

####################################

ani = animation.FuncAnimation(fig, updatefig, interval=500, blit=False, init_func=initfig)
plt.show()



