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

print "Plot density"

timestamp='003600'
nx=400
ny=400
xmin=0.
xmax=50.
ymin=0.
ymax=125.
filename='fields/' + timestamp +'.xy'

markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

# Read data, skipping 1st line
data=genfromtxt(filename,skip_header=1)
ng=data.shape[0]
numcols=data.shape[1]

print "ng:",ng,"numcols:",numcols

# extract
rhoe = data[:,0].reshape(ny,nx)
#rhoi = data[:,1].reshape(ny,nx)
#ex = data[:,2].reshape(ny,nx)
#ey = data[:,3].reshape(ny,nx)
#pot = data[:,4].reshape(ny,nx)

#print data
#print rhoi

fig = plt.figure(figsize=(8,8))
#fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.95) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels


#plt.subplot(121)
#cmap = plt.cm.PuOr # Set a colour map
#cmap.set_under ( 'w' ) # Low values set to 'w'hite
#cmap.set_bad ( 'w' ) # Bad values set to 'w'hite
plt.xlim((xmin, xmax))
plt.ylim((ymin, ymax))
im=plt.imshow(rhoi, cmap='PuOr_r',  interpolation='bicubic', origin='lower') 
plt.xlabel('$x$',fontsize=20)
plt.ylabel('$y$',fontsize=20)
plt.title('$n_e$',fontsize=20)
#cb = plt.colorbar(shrink=0.5)
#cb.set_label('rho e,i')


# contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)

	
plt.savefig(filename +'.eps') # Must occur before show()
plt.savefig(filename +'.png') # Must occur before show()
plt.show()
