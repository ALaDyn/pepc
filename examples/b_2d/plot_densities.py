#! /usr/bin/python

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




print "Plot field data"


plotboxsize   = 5.
animated = True
nx=100
ny=100


fig = plt.figure(figsize=(8,8))
fig.suptitle("Densities")
#fig.suptitle(filename, fontsize=26)
#fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels




def plot_image(field,position,color,ctitle):
    global fig

#    print field,position,color,ctitle
    plt.subplot(position)
    cmap = plt.get_cmap(color) # Set a colour map
    cmap.set_under ( 'w' ) # Low values set to 'w'hite
    cmap.set_bad ( 'w' ) # Bad values set to 'w'hite
    plt.imshow(field, cmap=cmap, aspect=1, interpolation='bilinear', vmin=0.01, vmax=1.5, origin='lower') 
    plt.xlabel('x')
    plt.ylabel('y')
    plt.colorbar()
#    cb.set_label(ctitle)
    plt.draw()
    return True

# contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)



def plot_from_file(fn,nx,ny):
    if os.path.isfile(fn):
        try:
            raw = []
            raw = genfromtxt(fn)
	    ng=raw.shape[0]
            numcols=raw.shape[1]
            print "ng:",ng,"numcols:",numcols
	# extract
	    rhoe = -raw[:,0].reshape(ny,nx)
            rhoi = raw[:,1].reshape(ny,nx)
	    plot_image(rhoe,211,'Reds','rhoe')
 	    plot_image(rhoi,212,'YlGn','rhoi')
            return
        except IOError:
	    print 'File ',fn,' not found'
            return False
    else:
	print 'File ',fn,' not found'
        return False


def plot_for_timestep(ts):
    global nx,ny
    filename = 'fields/%0*d'%(6, ts) + '.xy'
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

tmax=10000
increment = 100
plt.ion()
for timestamp in range(0,tmax,increment):
	plot_for_timestep(timestamp)
	sleep(0.1) # Time in seconds.
#	raw_input("Press key...")
	plt.clf()
#'	plt.show()
#	input = sys.stdin.readline() 

#plt.savefig(filename +'.png') # Must occur before show()
plt.show()
