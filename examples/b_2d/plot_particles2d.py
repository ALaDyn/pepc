#! /usr/bin/python

import matplotlib.pyplot as plt
#from matplotlib import rc
from numpy import *
#import matplotlib.font_manager
import matplotlib.cm as cm
import locale
import os

print "Plots phase space from particle dumps"

#locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
#rc('font', family='serif', size=28)
#rc('legend', fontsize='small')
#rc('xtick', labelsize='small')
#rc('ytick', labelsize='small')

filename='dumps/parts_p0001.000000'
FILE_COLUMN_X   =  0
FILE_COLUMN_Y   =  1
FILE_COLUMN_Z   =  2
FILE_COLUMN_Q   =  3
FILE_COLUMN_M   =  4
FILE_COLUMN_VX  =  5
FILE_COLUMN_VY  =  6
FILE_COLUMN_VZ  =  7
FILE_COLUMN_EX  =  8
FILE_COLUMN_EY  =  9
FILE_COLUMN_EZ  = 10
FILE_COLUMN_POT = 11
xmin=0
xmax=10
#markers=['+','*','.','1','2','3','4','<','>','D','H','^','_','d','h','o','p','s','v','x','|','None',' ','']
markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

data=genfromtxt(filename)
npart=data.shape[0]
numcols=data.shape[1]

print "npart:",npart,"numcols:",numcols
#print data
datanames=['x','y','z','vx','vy','vz','q','m','ex','ey','ez','pot','proc','label']
print "datanames:",datanames
fig = plt.figure(figsize=(8,8))
fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
#prop = matplotlib.font_manager.FontProperties(size=16)


####################################
ax=fig.add_subplot(111)
plt.xlabel(datanames[0])
plt.ylabel(datanames[1])
plt.grid(True, which='both')
col=data[:,6]/data[0,0]
print col[0:10]

ilim=npart/2
# electrons xy:
axScatter = plt.subplot(221)
axScatter.scatter(data[0:ilim,0],data[0:ilim,1],c='r',marker='o')
axScatter.set_aspect(1.)
axScatter.set_xlim(xmin,xmax)

# ions xy
axScatter = plt.subplot(222)
axScatter.scatter(data[ilim:npart,0],data[ilim:npart,1],c=col[ilim:npart],marker='o')
axScatter.set_aspect(1.)
axScatter.set_xlim(xmin,xmax)

# electrons x-px
axScatter = plt.subplot(223)
axScatter.scatter(data[0:ilim-1,0],data[0:ilim-1,3],c='r',marker='o')
#axScatter.set_aspect(1.)

# ions x-px
axScatter = plt.subplot(224)
axScatter.scatter(data[ilim:npart-1,0],data[ilim:npart-1,3],c=col[ilim:npart-1],marker='o')
#axScatter.set_aspect(1.)
#for colidx in range(1,7):
	#print colidx
#	plt.plot(data[:,0],data[:,colidx],label=datanames[colidx], linewidth=2.0,color=cm.hsv(24*colidx))#, markersize=5.0, markeredgewidth=1.0,marker=markers[colidx])
#	plt.legend(numpoints=1, ncol=1, bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

####################################

lims=plt.axis()
#plt.xticks(arange(lims[0],lims[1]+1,1.0)) 
#plt.yticks(0,1e3,200.)
#plt.minorticks_on()
	
plt.savefig(filename +'.png') # Must occur before show()
plt.show()
