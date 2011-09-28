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


currtimestep  = 1
stepincrement = 1
plotboxsize   = 5.
animated      = True


fig = plt.figure()
fig.suptitle("Radial Electrical Field")
ax = fig.gca()
maxy = 0.




def plot_particles_from_file(fn):
    if os.path.isfile(fn):
        try:
            raw = []
            raw = loadtxt(fn)
            return plot_radial_grid_fields(raw)
        except IOError:
            return False
    else:
        return False


def plot_radial_grid_fields(raw):
    global fig
    global maxy
    matplotlib.axes.Axes.cla(ax)

    npart     = len(raw[:,FILE_COLUMN_M])
    ngrid     = -1
    distance  = zeros(npart)
    field     = zeros(npart)
    potential = zeros(npart)
    for i in range(0,npart):
          if raw[i,FILE_COLUMN_M] == -1.:
              ngrid = ngrid + 1
              distance[ngrid] = sqrt(raw[i,FILE_COLUMN_X ]**2 + raw[i,FILE_COLUMN_Y ]**2 + raw[i,FILE_COLUMN_Z ]**2)
              field[ngrid]    = sqrt(raw[i,FILE_COLUMN_EX]**2 + raw[i,FILE_COLUMN_EY]**2 + raw[i,FILE_COLUMN_EZ]**2)
              potential[ngrid] = raw[i,FILE_COLUMN_POT]

    # prepare a line plot
    ax.plot(distance[1:ngrid-1], field[1:ngrid-1], linewidth=4)
    maxy = max(max(field[1:ngrid-1]),maxy)
    print maxy
    ax.set_xlim(0,plotboxsize)
    ax.set_xlabel('r')
    ax.set_ylabel('|E(r)|')
    ax.set_ylim(0.,maxy)
    fig.canvas.draw_idle()
    return True


def plot_particles_for_timestep(ts):
    filename = "snbs.output." + '%0*d'%(5, ts) + ".dat"
    if plot_particles_from_file(filename):
        print "Timestep: " + '%0*d'%(5, ts)
        ax.text(0,maxy, "Timestep: " + '%0*d'%(5, ts))
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
