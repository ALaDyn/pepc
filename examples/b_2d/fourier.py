#!/usr/bin/env python
"""
Example: simple Fourier analysis
Show how to make and save a simple line plot with labels, title and grid
"""
import numpy
from scipy import fftpack
import matplotlib.pyplot as plt

N=300  # sample rate
tmax=8.0 # 
dt = tmax/N #  sample interval
t = numpy.arange(0.0, tmax+dt, dt)  # time vector
theta = 2*numpy.pi*t  # norm. time (omega-t)
sigma=2.  # gaussian width
t0=4 # gaussian centre
s = (numpy.sin(5*theta)+0.25*numpy.cos(12*theta))*numpy.exp(-(t-t0)**2/sigma)


plt.subplot(211)
plt.plot(t, s)

plt.xlabel('time (s)')
plt.ylabel('voltage (mV)')
plt.title('Time series')
plt.grid(True)



n = N # length of the signal
k = numpy.arange(n)
wmax = numpy.pi/dt/2.
dw = 2*numpy.pi/tmax
wlim = 20.
frq = k*dw/2/numpy.pi # two sides frequency range
frq = frq[range(n/2)] # one side frequency range
print 'n,N,T,frq',n,N,dt
fs = fftpack.fft(s)/n # transform
fs = fs[range(n/2)] 
#w = numpy.arange(0.,wmax+dw, dw)

plt.subplot(212)
plt.yscale('log')
#ax.set_xscale('log')
plt.plot(frq, numpy.abs(fs), marker='o')
plt.xlim(0,wlim)
plt.ylim(1.e-8,1.)
plt.xlabel('Frequency (Hz)')
plt.ylabel('FT(U)')
plt.title('Frequency spectrum')
plt.grid(True)
plt.savefig('fourier')

plt.show()


