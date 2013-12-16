#! /usr/bin/python

import gobject
import matplotlib.pyplot as plt
#from matplotlib import rc
from numpy import *
import numpy as np
from scipy import fftpack
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


plotboxsize   = 3.
animated = True
tmax=5001
increment = 50
dt=0.1
nx=200
ny=200
xmin = 0.
xmax = 150.
ymin = 0.
ymax = 125.
xsample=35.
isamp = xsample/xmax*nx  # index of x-sample
dy = ymax/ny
y = dy*arange(1,ny+1,1)
#print y
fig = plt.figure(figsize=(6,6))


def plot_image(field,position,color,ctitle,fmin,fmax):
    global fig
    global xmin, xmax, ymin, ymax

#    print field,position,color,ctitle

    plt.subplot(position)
#    plt.subplots_adjust(bottom=0.1) 
    cmap = plt.get_cmap(color) # Set a colour map
#    cmap.set_under ( 'w' ) # Low values set to 'w'hite
#    cmap.set_bad ( 'w' ) # Bad values set to 'w'hite
    extent = [xmin, xmax, ymin, ymax]
    im=plt.imshow(field, cmap=cmap, aspect=1.6, extent=extent, interpolation='bilinear', origin='lower', vmin=fmin, vmax=fmax) 
    plt.xlabel('$x/\lambda_D$')
    plt.ylabel('$y/\lambda_D$')
    ax = im.get_axes()

#    plt.xticks(arange(xmin,xmax,25)) 
#    plt.yticks(arange(ymin,ymax,20))
    plt.minorticks_on()

    plt.colorbar(shrink=0.5)
    plt.title(ctitle,fontsize=15)
#    plt.set_label(ctitle)
    plt.savefig(filename +'.png')
    plt.draw()
    return True

# contour the gridded data, plotting dots at the nonuniform data points.
#CS = plt.contour(xi,yi,zi,15,linewidths=0.5,colors='k')
#CS = plt.contourf(xi,yi,zi,15,cmap=plt.cm.jet)



def plot_from_file(fn,nx,ny,is0,fmax):
    global y, dy, ymax, count, fgm
    if os.path.isfile(fn):
        try:
            data = []
            data = genfromtxt(fn,skiprows=1)
	    ng=data.shape[0]
            numcols=data.shape[1]
            print "ng:",ng,"numcols:",numcols
	# extract
	    rhoe = -data[:,0].reshape(ny,nx)
            rhoi = data[:,1].reshape(ny,nx)
	    ex = data[:,2].reshape(ny,nx)
	    ey = data[:,3].reshape(ny,nx)
	    pot = data[:,4].reshape(ny,nx)
	    plt.subplot(211)
	    plt.xlim(0,ymax)
    	    plt.grid(True)
	    specave=0.
	    nave=2

	    for isamp in range(is0-nave/2,is0+nave/2+1,1):
	    	f = ey[0:ny,isamp]


	    	plt.plot(y[1:ny],f[1:ny])
	    	dk = 2*np.pi/ymax
	    	k = dk*np.arange(ny)  # mode number range
	    	kmax = np.pi/dy
	    	klim = 20.
	    	k = k[range(ny/2)] # one side mode number range
#	    	print 'ny,kmax,',ny,kmax
	    	fs = fftpack.fft(f)/ny # transform
	    	fs = fs[range(ny/2)] 
	    	specave = specave + fs/(nave+1)


	    fmax = np.max(abs(specave))
  	    fgm[count]=fmax # store fastest growing mode
	    print fmax
	    plt.subplot(212)
	    plt.yscale('log')
#ax.set_xscale('log')
	    plt.plot(k, np.abs(specave))
	    plt.xlim(0,kmax)
	    plt.ylim(1.e-6,0.1)
	    plt.xlabel('Mode # ky')
	    plt.ylabel('$F(n_e)$')
	    plt.title('Spectrum')
	    plt.grid(True)


#	    plot_image(rhoe,221,'Reds','$n_e$',0.,1.0)
#  	    plot_image(pot,222,'BrBG','$\Phi$',-2.,2.)
# 	    plot_image(rhoi,222,'YlGn','$n_i$',0.,1.0)
# 	    plot_image(ex,223,'RdBu','Ex',-.5,.5)
# 	    plot_image(ey,224,'RdBu','Ey',-.5,.5)
            return True
        except IOError:
	    print 'File ',fn,' not found'
            return False
    else:
	print 'File ',fn,' not found'
        return False


def plot_for_timestep(ts):
    global nx,ny, filename, fig, isamp,fgm

    filename = 'fields/%0*d'%(6, ts) + '.xy'
    tsname = '$\omega_pt$=%0*d'%(6, ts)
    fig.suptitle(tsname)
    print filename,nx,ny
    if plot_from_file(filename,nx,ny,isamp,fmax):
        print "Timestep: " + '%0*d'%(6, ts)

        return True
    else:
        return False




#fn='fields/001000.xy'
#plot_from_file(fn,nx,ny,isamp)


plt.ion()
count=0
fgm = zeros((100))
for timestamp in range(500,2001,increment):
	count=count+1
	plot_for_timestep(timestamp)
#	sleep(0.1) # Time in seconds.
	raw_input("Press key...")
	plt.clf()


#'	plt.show()
#	input = sys.stdin.readline() 


#print fgm
t = dt*50*np.arange(count)
plt.plot(t[1:count], fgm[1:count], marker = 'o')
plt.xlim(0,dt*count)
plt.ylim(1.e-6,0.1)
plt.yscale('log')
plt.xlabel('t')
plt.ylabel('(kmax)')
plt.title('FGM')
plt.grid(True)
plt.savefig('ftey.png') # Must occur before show()
plt.show()
