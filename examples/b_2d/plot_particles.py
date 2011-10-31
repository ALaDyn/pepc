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


plotboxsize   = 10.
animated = True
nx=200
ny=100


fig = plt.figure(figsize=(10,10))
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



def plot_from_file(fn):
    global xmin,xmax,ymin,ymax
    if os.path.isfile(fn):
        try:
            data = []
            data = genfromtxt(fn)
	    nparts=data.shape[0]
            numcols=data.shape[1]
            print "np:",nparts,"numcols:",numcols
	    datanames=['x','y','z','vx','vy','vz','q','m','ex','ey','ez','pot','proc','label']
	    plt.xlabel(datanames[0])
            plt.ylabel(datanames[1])
            plt.grid(True, which='both')
	# extract
	    ilim=nparts/2
	    ilim=0
	# electrons xy:
 	    ax1 = plt.subplot(211)
	    ax1.scatter(data[0:ilim-1,0],data[0:ilim-1,1],c='r',marker='o')
	    ax1.scatter(data[ilim:nparts,0],data[ilim:nparts,1],c='b',marker='o')
	    ax1.set_aspect(1.)
	    ax1.set_xlim( (xmin, xmax) )
	    ax1.set_ylim( (ymin, ymax) )

 	    ax2 = plt.subplot(212)
	    ax2.scatter(data[0:ilim-1,0],data[0:ilim-1,3],c='r',marker='o')
	    ax2.scatter(data[ilim:nparts,0],data[ilim:nparts,3],c='b',marker='o')
	    ax2.set_aspect(1.)
	    ax2.set_xlim( (xmin, xmax) )
	    ax2.set_ylim( (-5.,5.) )
 	    plt.draw()
    	    return True
# ions xy
#axScatter = plt.subplot(222)
#axScatter.scatter(data[ilim:npart-1,0],data[ilim:npart-1,1],c=col[ilim:npart-1],marker='o')
#axScatter.set_aspect(1.)

        except IOError:
	    print 'File ',fn,' not found'
            return False
    else:
	print 'File ',fn,' not found'
        return False


def plot_for_timestep(ts):

    filename = 'dumps/parts_p0000.%0*d'%(6, ts)
    if plot_from_file(filename):
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

tmax=100
increment = 1
xmin=-10
xmax=20
ymin=0
ymax=10

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
