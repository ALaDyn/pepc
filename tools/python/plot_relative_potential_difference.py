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


fig = plt.figure()
#fig.suptitle("Simulated Particles")
ax = fig.gca()
ax.set_xlabel('r')
ax.set_ylabel('relative error of potential')


rawpepc = []
rawpepc = loadtxt("pepc.output")
rawdirect = []
rawdirect = loadtxt("direct.output")

dist = zeros(len(rawpepc[:,FILE_COLUMN_POT]))
dp   = zeros(len(rawpepc[:,FILE_COLUMN_POT]))

diffpot = abs(rawpepc[:,FILE_COLUMN_POT] - rawdirect[:,FILE_COLUMN_POT])/abs(rawdirect[:,FILE_COLUMN_POT])

i=-1
for p in range(0,len(rawpepc[:,FILE_COLUMN_X])):
    if rawpepc[p,FILE_COLUMN_M] == -2.:
        i = i+1
        dist[i] = sqrt(rawpepc[p,FILE_COLUMN_X]**2 + rawpepc[p,FILE_COLUMN_Y]**2 + rawpepc[p,FILE_COLUMN_Z]**2)
        dp[i]   = diffpot[p]

ax.semilogy(dist[0:i], dp[0:i])

plt.show()
