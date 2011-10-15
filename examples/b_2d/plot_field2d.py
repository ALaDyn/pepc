#! /usr/bin/python

import matplotlib.pyplot as plt
#from matplotlib import rc
from numpy import *
import numpy as np
#import matplotlib.font_manager
from matplotlib.mlab import griddata
import matplotlib.delaunay
import matplotlib.cm as cm
import locale
import os

print "Plot field data"

timestamp='000150'
nx=200
ny=100
filename='fields/' + timestamp +'.xy'

markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

data=genfromtxt(filename)
ng=data.shape[0]
numcols=data.shape[1]

print "ng:",ng,"numcols:",numcols

# extract
rhoe = -data[:,0].reshape(ny,nx)
rhoi = data[:,1].reshape(ny,nx)

#print data
#print rhoi

fig = plt.figure(figsize=(8,8))
#fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels


plt.subplot(211)
cmap = plt.cm.Reds # Set a colour map
cmap.set_under ( 'w' ) # Low values set to 'w'hite
cmap.set_bad ( 'w' ) # Bad values set to 'w'hite
plt.imshow(rhoe, cmap=cmap, aspect=1, interpolation='bilinear', vmin=0.01, origin='lower') 
plt.xlabel('x')
plt.ylabel('y')
cb = plt.colorbar()
cb.set_label('rhoe')

plt.subplot(212)
cmap = plt.cm.YlGn # Set a colour map
cmap.set_under ( 'w' ) # Low values set to 'w'hite
cmap.set_bad ( 'w' ) # Bad values set to 'w'hite
plt.imshow(rhoi, cmap=cmap, aspect=1, interpolation='bilinear', vmin=0.01, origin='lower') 
plt.xlabel('x')
plt.ylabel('y')

#plt.axis([x.min(), x.max(), y.min(), y.max()])
cb = plt.colorbar()
cb.set_label('rhoi')

# contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)

	
plt.savefig(filename +'.png') # Must occur before show()
plt.show()
