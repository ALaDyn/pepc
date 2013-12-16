#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc


filename = 'energy.dat'
columns = [ '$t / {\\rm fs}$',
            '$U_{\\rm total}^{\\rm (pot)}$',
            '$U_{\\rm near field}^{\\rm (pot)}$',
            '$U_{\\rm far field}^{\\rm (pot)}$',
            '$U_{\\rm e}^{\\rm (kin)}$',
            '$U_{\\rm i}^{\\rm (kin)}$',
            '$U_{\\rm tot}^{\\rm (kin)}$',
            '$U_{\\rm e}^{\\rm (kin, wo drift)}$',
            '$U_{\\rm i}^{\\rm (kin, wo drift)}$',
            '$U^{\\rm (tot)}$',
            '$T_{\\rm e}^{\\rm (0)}$',
            '$T_{\\rm e}^{\\rm (uncorr)}$',
            '$\chi_{\\rm e}$',
            '$\delta T_{\\rm i}$',
            '$T_{\\rm i}^{\\rm (0)}$',
            '$T_{\\rm i}^{\\rm (uncorr)}$',
            '$\chi_{\\rm i}$',
            '$\delta T_{\\rm i}$' ]
ylab = '$E / {\\rm Ryd}$'


data = np.loadtxt(filename)

rc('text', usetex=False)
rc('font', family='serif', size=12)
fig = plt.figure(1, figsize=(10,10))

plt.suptitle(filename, fontsize=30)

ax1 = plt.subplot(211)
plt.plot(data[:,0],data[:,4],data[:,0],data[:,5],data[:,0],data[:,6])
ax1.set_xlabel(columns[0], fontsize=20)
ax1.set_ylabel(ylab, fontsize=20)
leg = ax1.legend((columns[4],columns[5],columns[6]),
           'lower right', shadow=False, frameon=False,ncol=3)

ax2 = plt.subplot(212)
plt.plot(data[:,0],data[:,1],data[:,0],data[:,9])
ax1.set_xlabel(columns[0], fontsize=20)
ax1.set_ylabel(ylab, fontsize=20)
leg = ax2.legend((columns[1],columns[9]),
           'lower right', shadow=False, frameon=False,ncol=3)


plt.savefig('energy.pdf')
plt.show()

