#! /usr/bin/python

import gobject
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
from numpy import *
import numpy as np
#import matplotlib.font_manager
from matplotlib.mlab import griddata
import matplotlib.delaunay
import matplotlib.cm as cm
from matplotlib.colors import LogNorm
import locale
import os
import os.path
import sys
from time import sleep
#rc('text',usetex=True)
#rcParams.update({'font.size': 16})




print "Plot field data"


plotboxsize   = 5.
animated = True
tmax=410
increment = 20
dt=0.15
nx=200
ny=200
xmin = 0.
xmax = 100.
ymin = 0.
ymax = 125.

fig = plt.figure(figsize=(10,10))
#fig.suptitle("Fields")
#fig.suptitle(filename, fontsize=26)
#fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels




def plot_image(field,position,color,ctitle,fmin,fmax):
    global fig
    global xmin, xmax, ymin, ymax

#    print field,position,color,ctitle

    plt.subplot(position)
#    plt.subplots_adjust(bottom=0.1) 
    cmap = plt.get_cmap(color) # Set a colour map
    lowest = cmap(fmin) # returns RGBA 4-tuple
#    print 'map bottom:',lowest
    cmap.set_under ( lowest ) # Low values set to cmap min
    cmap.set_bad ( lowest ) # Bad values set to cmap min
    extent = [xmin, xmax, ymin, ymax]
#  linear scale
#    im=plt.imshow(field, cmap=cmap, aspect=1., extent=extent, interpolation='bilinear', origin='lower', vmin=fmin, vmax=fmax) 
#  log scale
    im=plt.imshow(field, cmap=cmap, aspect=1., extent=extent, interpolation='bilinear', origin='lower', norm=LogNorm(vmin=fmin, vmax=fmax))
    plt.xlabel(r'$x/\lambda_D$')
    plt.ylabel(r'$y/\lambda_D$')
    ax = im.get_axes()

#    plt.xticks(arange(xmin,xmax,25)) 
#    plt.yticks(arange(ymin,ymax,20))
    plt.minorticks_on()

#    plt.colorbar(shrink=0.5)
    plt.title(ctitle,fontsize=15)
#    plt.set_label(ctitle)
    plt.savefig(filename +'.png')
    plt.savefig(filename +'.pdf')
#    plt.draw()
    return True

# contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)



def plot_from_file(fn,nx,ny):
    if os.path.isfile(fn):
        try:
            data = []
            data = genfromtxt(fn,skip_header=1)
	    ng=data.shape[0]
            numcols=data.shape[1]
#            print "ng:",ng,"numcols:",numcols
	# extract
	    rhoe = -data[:,0].reshape(ny,nx)
#            rhoi = data[:,1].reshape(ny,nx)
#	    ex = data[:,2].reshape(ny,nx)
#	    ey = data[:,3].reshape(ny,nx)
#	    pot = data[:,4].reshape(ny,nx)
	    plot_image(rhoe,111,'PuOr_r','r$n_e$',0.0001,1.5)
#  	    plot_image(pot,222,'BrBG','$\Phi$',-2.,2.)
# 	    plot_image(rhoi,222,'YlGn','$n_i$',0.,1.0)
# 	    plot_image(ex,223,'RdBu','Ex',-.2,.2)
# 	    plot_image(ey,224,'RdBu','Ey',-.2,.2)
#            return
        except IOError:
	    print 'File ',fn,' not found'
            return False
    else:
	print 'File ',fn,' not found'
        return False


def plot_for_timestep(ts):
    global nx,ny, filename, fig, dt

    filename = 'fields/%0*d'%(6, ts) + '.xy'
    tsname = '$\omega_pt$=%0*d'%(6, ts*dt)
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


# Execute update method every 500ms
#if animated:
    #gobject.idle_add(next_plot)
#    gobject.timeout_add(250, next_plot)
#else:
#    next_plot()

#gobject.idle_add(next_plot)

#fn='fields/000000.xy'
#plot_from_file(fn,nx,ny)


plt.ion()
for timestamp in range(0,tmax,increment):
	plot_for_timestep(timestamp)
#	sleep(0.1) # Time in seconds.
#	raw_input("Press key...")
	plt.clf()
#'	plt.show()
#	input = sys.stdin.readline() 

#plt.savefig(filename +'.png') # Must occur before show()
#plt.show()
