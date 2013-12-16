#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import *
import matplotlib.font_manager
import matplotlib.cm as cm
import locale
import os.path
from pylab import *


rc('text', usetex=True)
rc('font', family='serif', size=28)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

FIELD_NPARTS=3
FIELD_PE=1
FIELD_FETCHES=9
FIELD_SHIPS=10
FIELD_INTERACTIONS=11
FIELD_MACEVALS=12
FIELD_RELWORK=13

fileprefix='stats.'

fieldnames=["Number of Particles", "Number of Fetches", "Number of Ships", "Number of Interactions", "Number of Mac Evaluations", "relative Work"]
fields=[3, 9, 10, 11, 12, 13]

for i in range(0,len(fields)):
	field=fields[i]
	filein='extr.stats.FIELD' + `field` + '.dat'
	
	raw=[]
	raw=transpose(loadtxt(filein))
	
	fig = plt.figure(figsize=(12,8))
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
	ax.set_title(fieldnames[i])
	plt.xlabel('Timestep')
	plt.ylabel('MPI rank')
	
	im = ax.imshow(raw)
	fig.colorbar(im)
	plt.savefig(fileprefix + fieldnames[i].replace(' ', '') + '.pdf')
	#plt.show()
