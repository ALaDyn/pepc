#! /usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import *
import matplotlib.font_manager
import matplotlib.cm as cm
import locale
import os

print "Plots data from energy.dat"
print "optimized for pepc-mw frontend"
print "Creating energy.dat.pdf in current directory"

locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

rc('text', usetex=True)
rc('font', family='serif', size=28)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

filename='energy.dat'

#markers=['+','*','.','1','2','3','4','<','>','D','H','^','_','d','h','o','p','s','v','x','|','None',' ','']
markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

data=loadtxt(filename)
numlines=data.shape[0]
numcols=data.shape[1]

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

fig = plt.figure(figsize=(12,8))
fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
prop = matplotlib.font_manager.FontProperties(size=16)


####################################
ax=fig.add_subplot(111)
plt.xlabel(datanames[0])
plt.ylabel("Energy [Ry]")
plt.grid(True, which='both')
#ax.set_yscale('log')

for colidx in range(1,10):
	plt.plot(data[:,0],data[:,colidx],label=datanames[colidx], linewidth=1.0,color=cm.hsv(24*colidx))#, markersize=1.0, markeredgewidth=1.0,marker=markers[colidx])
plt.legend(numpoints=1, prop=prop, ncol=1, frameon=True, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

####################################

lims=plt.axis()
#plt.xticks(arange(lims[0],lims[1],2.0)) 
plt.minorticks_on()
	
plt.savefig(filename +'.pdf') # Must occure before show()
plt.show()
