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




print "Particle tracks"
plotboxsize   = 8.
animated = True


fig = plt.figure(figsize=(10,10))
fig.suptitle("Tracks")
#fig.suptitle(filename, fontsize=26)
#fig.subplots_adjust(right=0.8) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels



def plot_from_file(fn):
    global xmin,xmax,ymin,ymax
    global ax1
    print ax1
    if os.path.isfile(fn):
        try:
            data = []
            data = genfromtxt(fn)
	    nparts=data.shape[0]
            numcols=data.shape[1]
            print "np:",nparts,"numcols:",numcols
	# extract
	    ilim=nparts/2
	# electrons xy:
	    ax1.scatter(data[0:ilim-1,0],data[0:ilim-1,1],c='r',marker='o')
	    ax1.scatter(data[ilim:nparts,0],data[ilim:nparts,1],c='b',marker='o')
 	    plt.draw()
    	    return True

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

tmax=101
increment = 1
xmin=0
xmax=10
ymin=0
ymax=10

plt.ion()

datanames=['x','y','z','vx','vy','vz','q','m','ex','ey','ez','pot','proc','label']
plt.xlabel(datanames[0])
plt.ylabel(datanames[1])
plt.grid(True, which='both')
ax1 = plt.subplot(111)
ax1.set_xlim( (xmin, xmax) )
ax1.set_ylim( (ymin, ymax) )
ax1.set_aspect(1.)
plt.draw()

for timestamp in range(0,tmax,increment):
#	plt.clf()
	plot_for_timestep(timestamp)
#	sleep(0.1) # Time in seconds.
#	raw_input("Press key...")
#'	plt.show()
#	input = sys.stdin.readline() 

#plt.savefig(filename +'.png') # Must occur before show()
plt.show()
