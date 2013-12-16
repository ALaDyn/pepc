#!/usr/bin/python

import gobject
from numpy import *
from mpl_toolkits.mplot3d import Axes3D
import matplotlib
matplotlib.use('GTKAgg') # do this before importing pylab
import matplotlib.pyplot as plt
import os.path


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


currtimestep  = 0
stepincrement = 1
plotboxsize   = 5.
animated      = True


fig = plt.figure()
fig.suptitle("Simulated Particles")
ax = fig.gca(projection='3d')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.axis('equal')
cm = matplotlib.cm.get_cmap('jet')



def plot_particles_from_file(fn):
    if os.path.isfile(fn):
        try:
            raw = []
            raw = loadtxt(fn)
            return plot_raw_data(raw)
        except IOError:
            return False
    else:
        return False



def plot_raw_data(raw):
    global fig
    matplotlib.axes.Axes.cla(ax)
    ax.scatter(raw[:,FILE_COLUMN_X],  raw[:,FILE_COLUMN_Y],  raw[:,FILE_COLUMN_Z], cmap=cm, c=raw[:,FILE_COLUMN_Q], s=5, linewidth=0, vmax=0.01, vmin=-0.01)
    plt.axes().set_aspect('equal', 'datalim')
    ax.set_xlim3d(-plotboxsize,plotboxsize)
    ax.set_ylim3d(-plotboxsize,plotboxsize)
    ax.set_zlim3d(-plotboxsize,plotboxsize)
    #ax.mouse_init()
    fig.canvas.draw()
    return True


def plot_particles_for_timestep(ts):
    filename = "parts_p0000." + '%0*d'%(5, ts)
    if plot_particles_from_file(filename):
        print "Timestep: " + '%0*d'%(5, ts)
        return True
    else:
        return False



def next_plot():
    global currtimestep

    if plot_particles_for_timestep(currtimestep):
        currtimestep = currtimestep + stepincrement

    return True


# Execute update method every 500ms
if animated:
    #gobject.idle_add(next_plot)
    gobject.timeout_add(250, next_plot)
else:
    next_plot()

# Display the plot
plt.show()
