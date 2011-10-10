#! /usr/bin/python

import matplotlib.pyplot as plt
#from matplotlib import rc
from numpy import *
#import matplotlib.font_manager
import matplotlib.cm as cm
import locale
import os

print "Plots data from simple.dat"
#print "Creating energy.dat.pdf in current directory"

#locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
#rc('font', family='serif', size=28)
#rc('legend', fontsize='small')
#rc('xtick', labelsize='small')
#rc('ytick', labelsize='small')

filename='energy.dat'

#markers=['+','*','.','1','2','3','4','<','>','D','H','^','_','d','h','o','p','s','v','x','|','None',' ','']
markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

#data=genfromtxt(filename)
data=loadtxt(filename,skiprows=1)
numlines=data.shape[0]
numcols=data.shape[1]

print "numlines:",numlines,"numcols:",numcols
#print data
datanames=['t','Upot','Umag','Ue','Ui','Ub','Utot']
print "datanames:",datanames
fig = plt.figure(figsize=(12,8))
fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
#prop = matplotlib.font_manager.FontProperties(size=16)


####################################
ax=fig.add_subplot(111)
plt.xlabel(datanames[0])
plt.grid(True, which='both')
#ax.set_yscale('log')

#plt.plot(data[:,0],data[:,1],label=datanames[1], linewidth=1.0,color=cm.hsv(24*1), markersize=5.0, markeredgewidth=1.0,marker=markers[1])

for colidx in range(1,7):
	#print colidx
	plt.plot(data[:,0],data[:,colidx],label=datanames[colidx], linewidth=2.0,color=cm.hsv(24*colidx))#, markersize=5.0, markeredgewidth=1.0,marker=markers[colidx])
	plt.legend(numpoints=1, ncol=1, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

####################################

lims=plt.axis()
plt.xticks(arange(lims[0],lims[1]+1,1.0)) 
#plt.yticks(0,1e3,200.)
plt.minorticks_on()
	
plt.savefig(filename +'.pdf') # Must occur before show()
plt.show()
