#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


filename = 'cluster.dat'
columns = ['itime', '$t / {\\rm fs}$', '$r_{\\rm cluster}^{\\rm (rms)} / a_{\\rm B}$', '$N_{\\rm ion}$', '$N_{\\rm e}^{\\rm (free)}$', '$N_{\\rm e}^{\\rm (bound)}$', '$N_{\\rm  e}^{(r<r_{\\rm cluster}^{\\rm (rms)})}$', '$N_{\\rm e}^{\\rm (E<0)}$']


data = np.loadtxt(filename)

rc('text', usetex=False)
rc('font', family='serif', size=12)
fig = plt.figure(1, figsize=(10,10))

plt.suptitle(filename, fontsize=30)

ax1 = plt.subplot(211)
plt.plot(data[:,1],data[:,2])
ax1.set_xlabel(columns[1], fontsize=20)
ax1.set_ylabel(columns[2], fontsize=20)

ax2 = plt.subplot(212)
plt.plot(data[:,1],data[:,5],'k',data[:,1],data[:,6],'g',data[:,1],data[:,7],'r')
ax2.set_xlabel(columns[1], fontsize=20)
ax2.set_ylabel(columns[5], fontsize=20)
leg = ax2.legend((columns[5],columns[6],columns[7]),
           'upper right', shadow=False, frameon=False,ncol=3)


plt.savefig('cluster.pdf')
plt.show()

