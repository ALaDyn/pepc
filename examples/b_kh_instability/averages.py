#! /usr/bin/python

import gobject
import matplotlib.pyplot as plt
#from matplotlib import rc
from numpy import *
import numpy as np
#from scipy import fftpack
#import matplotlib.font_manager
from matplotlib.mlab import griddata
import matplotlib.delaunay
import matplotlib.cm as cm
import locale
import os
import os.path
import sys
from time import sleep

print "Extract y-averaged lineouts"

directory='fields'
plotboxsize   = 3.
animated = True
tmin=0
tmax=401
increment = 20
dt=0.15
nx=200
ny=200
xmin = 0.
xmax = 100.
ymin = 0.
ymax = 125.
nmax = 1.5
dx = xmax/nx
x = dx*arange(1,nx+1,1) # x-axis
dy = ymax/ny
y = dy*arange(1,ny+1,1)
#print y
fig = plt.figure(figsize=(8,8))


#############################




def plot_from_file(fn,nx,ny,is0,fmax):
    global y, dy, ymax, x,dx,xmax, count, fgm, ne, ni, exav
    if os.path.isfile(fn):
        try:
            data = []
            data = genfromtxt(fn,skip_header=1)
	    ng=data.shape[0]
            numcols=data.shape[1]
#            print "ng:",ng,"numcols:",numcols
	# extract
	    rhoe = -data[:,0].reshape(ny,nx)
            rhoi = data[:,1].reshape(ny,nx)
	    ex = data[:,2].reshape(ny,nx)
#	    ey = data[:,3].reshape(ny,nx)
#	    pot = data[:,4].reshape(ny,nx)

      	    ne = rhoe.mean(0) # Average over y-axis
      	    ni = rhoi.mean(0) # Average over y-axis
      	    exav = ex.mean(0) # Average over y-axis

	    plt.subplot(311)
#	    plt.yscale('log')
	    plt.plot(x, ne)
	    plt.plot(x, ni)
	    plt.xlim(0,xmax) 
	    plt.ylim(0,fmax)
#	    plt.xlabel('$x/\lambda_D$',fontsize=15)
	    plt.ylabel('$<n_{e,i}>$',fontsize=15)
	    plt.grid(True)

	    plt.subplot(312)
	    plt.yscale('log')
	    plt.plot(x, ne)
	    plt.plot(x, ni)
	    plt.xlim(0,xmax) 
	    plt.ylim(1.e-5,fmax)
#	    plt.xlabel('$x/\lambda_D$',fontsize=15)
	    plt.ylabel('$<n_{e,i}>$',fontsize=15)
	    plt.grid(True)

	    plt.subplot(313)
#	    plt.yscale('log')
	    plt.plot(x, exav)
	    plt.xlim(0,xmax) 
	    plt.ylim(-0.8,0.8)
	    plt.xlabel('$x/\lambda_D$',fontsize=15)
	    plt.ylabel('$<Ex>$',fontsize=15)
	    plt.grid(True)


            return True
        except IOError:
	    print 'File ',fn,' not found'
            return False
    else:
	print 'File ',fn,' not found'
        return False


class prettyfloat(float):
    def __repr__(self):
        return "%0.2f" % self

def plot_for_timestep(ts):
    global nx,ny, directory, fig, isamp,fgm,nmax, dt,x,ne,ni,exav

    filename = directory+'/%0*d'%(6, ts) + '.xy'
    tsname = '%0*d'%(6, ts*dt)
    wpts = '$\omega_pt$='+tsname
    fig.suptitle(directory+':  '+wpts)
#    print filename,nx,ny
    if plot_from_file(filename,nx,ny,isamp,nmax):
        print "Timestep: " + '%0*d'%(6, ts), filename,nx,ny
# Write png file
	plt.savefig(directory+'/fave.'+tsname+'.png') # Must occur before show()

# Write out y-averages to ascii file:
        FILE = open(directory+'/fave.'+tsname,"w")

	for ix in range(1,nx,1):
	      cx="{0:0.2f}".format(x[ix])
	      cne="{0:0.2e}".format(ne[ix])
	      cni="{0:0.2e}".format(ni[ix])
	      cexav="{0:0.2e}".format(exav[ix])
	      cline=str(cx)+"\t"+str(cne).rjust(10)+"\t"+str(cni).rjust(10)+"\t"+str(cexav).rjust(10)+"\n"
    	      FILE.write(cline)
    
        FILE.close()


        return True
    else:
        return False


def plot_series(fs,ctitle):
	t = dt*50*np.arange(count)
	plt.plot(t[1:count], fs[1:count], marker = 'o')
	plt.xlim(0,dt*count)
	plt.ylim(1.e-6,0.1)
	plt.yscale('log')
	plt.xlabel('t')
	plt.ylabel('(kmax)')
	plt.title(ctitle)
	plt.grid(True)
	plt.savefig('ftey.png') # Must occur before show()



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

########################

fn='fields7/000100.xy'
isamp=1
#plot_from_file(fn,nx,ny,isamp,nmax)
#plt.show()

plt.ion()
count=0

for timestamp in range(tmin,tmax,increment):
	count=count+1
	plot_for_timestep(timestamp)
#	plt.draw()
#	sleep(0.1) # Time in seconds.
#	raw_input("Press key...")

	plt.clf()

#plt.show()
#	    plot_image(rhoe,221,'Reds','$n_e$',0.,1.0)
#plot_series(fgm)

