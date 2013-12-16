#!/usr/bin/env python
"""
Example: simple line plot.
Show how to make and save a simple line plot with labels, title and grid
"""
import numpy
import pylab
a=1.0
eps=0.8*a  # norm to a
np=100
inc=4./np
r = numpy.arange(0.0, 4.0+inc, inc)
d = r
ieps = int(eps/inc)
f=numpy.ones(np)
for i in range(ieps,np,1):
	f[i] =2.*a**8/d[i]**8 - 1.*a**4/d[i]**4

ymax = numpy.max(f)
f[1:ieps]=ymax
print ieps,ymax
#print f
pylab.plot(d[1:np], f[1:np])

pylab.xlabel('r ')
pylab.ylabel('force ')
pylab.xlim((0.,6.))
pylab.ylim((-1.,2*ymax))
pylab.title('Lennard-Jones')
pylab.grid(True)
pylab.savefig('lj_plot')

pylab.show()


