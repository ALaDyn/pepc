#! /usr/bin/python
# 
# Script to analyse Kelvin-Helmholtz instability developing
# in magnetised plasma sheath.
# Takes FT of field quantities in y, averaged over several strips in x-direction

import gobject
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
import os.path
import sys
from time import sleep

print "Fourier analysis of fields"

plotboxsize   = 5.
animated = True
tmax=51
increment = 50
nx=200     # No. grid points in 2D plots
ny=200
xmin = 0.  # Actual box dimensions
xmax = 150.
ymin = 0.
ymax = 100.

fig = plt.figure(figsize=(10,10))
#fig.suptitle("Fields")
#fig.suptitle(filename, fontsize=26)



def plot_image(field,position,color,ctitle,fmin,fmax):
    global fig
    global xmin, xmax, ymin, ymax

#    print field,position,color,ctitle

    plt.subplot(position)
#    plt.subplots_adjust(bottom=0.1) 
    cmap = plt.get_cmap(color) # Set a colour map
#    cmap.set_under ( 'w' ) # Low values set to 'w'hite
#    cmap.set_bad ( 'w' ) # Bad values set to 'w'hite
    extent = [xmin, xmax, ymin, ymax]
    im=plt.imshow(field, cmap=cmap, aspect=1.6, extent=extent, interpolation='bilinear', origin='lower', vmin=fmin, vmax=fmax) 
    plt.xlabel('$x/\lambda_D$')
    plt.ylabel('$y/\lambda_D$')
    ax = im.get_axes()

#    plt.xticks(arange(xmin,xmax,25)) 
#    plt.yticks(arange(ymin,ymax,20))
    plt.minorticks_on()

    plt.colorbar(shrink=0.5)
    plt.title(ctitle,fontsize=15)
#    plt.set_label(ctitle)
    plt.savefig(filename +'.png')
    plt.draw()
    return True

# contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)



def plot_from_file(fn,nx,ny):
    if os.path.isfile(fn):
        try:
            data = []
            data = genfromtxt(fn,skiprows=1)
	    ng=data.shape[0]
            numcols=data.shape[1]
            print "ng:",ng,"numcols:",numcols
	# extract
	    rhoe = -data[:,0].reshape(ny,nx)
            rhoi = data[:,1].reshape(ny,nx)
	    ex = data[:,2].reshape(ny,nx)
	    ey = data[:,3].reshape(ny,nx)
	    pot = data[:,4].reshape(ny,nx)
	    plot_image(rhoe,221,'Reds','$n_e$',0.,1.0)
 	    plot_image(rhoi,222,'YlGn','$n_i$',0.,1.0)
 	    plot_image(ex,223,'RdBu','Ex',-.5,.5)
 	    plot_image(ey,224,'RdBu','Ey',-.5,.5)
#  	    plot_image(pot,325,'BrBG','$\Phi$',-5.,5.)
            return
        except IOError:
	    print 'File ',fn,' not found'
            return False
    else:
	print 'File ',fn,' not found'
        return False


def plot_for_timestep(ts):
    global nx,ny, filename, fig

    filename = 'fields/%0*d'%(6, ts) + '.xy'
    tsname = '$\omega_pt$=%0*d'%(6, ts)
    fig.suptitle(tsname)
    print filename,nx,ny
    if plot_from_file(filename,nx,ny):
        print "Timestep: " + '%0*d'%(6, ts)
        return True
    else:
        return False



def next_plot():
    global timestamp
    print timestamp
    if plot_for_timestep(timestamp):
        timestamp = timestamp + increment

    return True


plt.ion()
for timestamp in range(0,tmax,increment):
	plot_for_timestep(timestamp)
#	sleep(0.1) # Time in seconds.
#	raw_input("Press key...")
	plt.clf()
#'	plt.show()
#	input = sys.stdin.readline() 

#plt.savefig(filename +'.png') # Must occur before show()
plt.show()