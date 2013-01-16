#! /lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import *
import matplotlib.font_manager
import matplotlib.cm as cm
import locale
import os
import sys
import json
from scipy import polyval, polyfit

locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
rc('font',   family='serif', size=24)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

filename='nose_hoover.dat'

includefiles = False
if len(sys.argv) > 1:
  subdir = sys.argv[1]
  includefiles = True

print "Usage: %s [datadir with PRE*.dat-files]" % sys.argv[0]

#markers=['+','*','.','1','2','3','4','<','>','D','H','^','_','d','h','o','p','s','v','x','|','None',' ','']
markers=['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

data=loadtxt(filename)
numlines=data.shape[0]
numcols=data.shape[1]

datanames=['$t [\mathrm{fs}]$',
           '$U_{\mathrm{pot}}^{\mathrm{(e)}}$',
           '$T_{\mathrm{new}}^{\mathrm{(e)}}$',
           '$T_{\mathrm{0}}^{\mathrm{(e)}}$',
           '$\eta^{\mathrm{(e)}}$',
           '$v_{\eta^{\mathrm{(e)}}}$',
           '$Q^{\mathrm{(e)}}$',
           '$\mathcal{H}^{\mathrm{(e)}}$',
           '$U_{\mathrm{pot}}^{\mathrm{(i)}}$',
           '$T_{\mathrm{new}}^{\mathrm{(i)}}$',
           '$T_{\mathrm{0}}^{\mathrm{(i)}}$',
           '$\eta^{\mathrm{(i)}}$',
           '$v_{\eta^{\mathrm{(i)}}}$',
           '$Q^{\mathrm{(i)}}$',
           '$\mathcal{H}^{\mathrm{(i)}}$',
           '$\mathcal{H}^{\mathrm{(e)}}+\mathcal{H}^{\mathrm{(i)}}$']

fig = plt.figure(figsize=(16,12))
fig.suptitle(filename, fontsize=20)
fig.subplots_adjust(right=0.95, top=0.9, bottom=0.1,left=0.1) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
prop = matplotlib.font_manager.FontProperties(size=12)


####################################
if os.path.exists("./simparams.dat"):
  f = open('./simparams.dat', 'r')
  simparams = json.load(f)
  f.close()
else:
  simparams = {'defaultfile' : True, 't0': 3.5, 'w_over_wpl' : 3.0, 'E0' : 1.56542, 'wpl' : 0.272921, 'Z' : 1, 'ne' : 0.00148185, 'v0' : 5.40775, 'Ne' : 500}
  f = open('./simparams.dat', 'w')
  json.dump(simparams, f)
  f.close()

if simparams['defaultfile']:
  print ">>>>>>>>>> simparams.dat contains default vaules - you should correct this and set defaultfile=False therein before continuing <<<<<<<<<<"
  sys.exit(0)
  
# we skip everything before switching on the laser
for startidx in range(len(data[:,0])):
  if data[startidx,0] > simparams['t0']:
    break
    
    
###################################
# linear fit for He+Hi  wit   y = a*x + b
FITCOL = 15
(a, b) = polyfit(data[startidx:-1,0], data[startidx:-1, FITCOL], 1)

print "Fit result: %f * t + %f" % (a, b)

yfit   = polyval( [a,b], data[startidx:-1,0])


###################################
ax = fig.add_subplot(2,2,1)
ax.set_xlabel(datanames[0])
ax.set_ylabel(r'$\mathcal{H}\,\mathrm{[Ry]}$')
ax.grid(True, which='both')
#ax.set_yscale('log')
#ax.set_xlim([0.,30.])

for p in [7, 14, 15]:
  ax.plot(data[startidx:-1,0], data[startidx:-1,p],  '-', label=datanames[p], linewidth=1.5)

ax.plot(data[startidx:-1,0], yfit,  '--', label="%s - fit" % datanames[FITCOL], linewidth=1.0)

ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='upper left', bbox_to_anchor=(1.02, 0.0, 1.3, 1.0))
plt.minorticks_on()

###################################
ax = fig.add_subplot(2,2,3)
ax.set_xlabel(datanames[0])
ax.set_ylabel("Temperature [Ry]")
ax.grid(True, which='both')
#ax.set_yscale('log')
#ax.set_xlim([0.,30.])

for p in [2, 3, 9, 10]:
  ax.plot(data[startidx:-1,0], data[startidx:-1,p],  '-', label=datanames[p], linewidth=1.5)

ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='right')


####################################
ax = fig.add_subplot(2,2,4)

kB   = 1.
eps0 = 1./(4.*pi)
unit_t0_in_s       =       4.8377687E-17
unit_t0_in_fs      = unit_t0_in_s / 1.E-15
unit_qe            = -sqrt(2.) 

# function for temp-dependent collision frequency
def nu(Temp, dU_dT, s):
        if 'v0_over_vtherm' in s.keys():
          # Use expression for fixed v0/vtherm
          Etherm = 3./2.*kB*Temp
          return 1./s['v0_over_vtherm'] * 4. / Etherm * dU_dT / s['wpl']
        else:
          # Use expression for fixed E0
          return s['w_over_wpl']**2 * 2./(eps0*s['E0']**2) * s['ne'] * dU_dT / s['wpl']


Tempe = data[-1,2]

ax.set_xlabel(r'$\Gamma$')
ax.set_ylabel(r'$\nu_{\rm ei} / \omega_{\rm pl}$')
ax.grid(True, which='both')

Gamma = simparams['Z']*unit_qe**2/(kB*Tempe)*(4.*pi/3.*simparams['ne'])**(1./3.)
nuval = nu(Tempe, a*unit_t0_in_fs/(simparams["Ne"] + simparams["Ni"]), simparams)

ax.set_xscale('log')
ax.set_yscale('log')

ax.plot(Gamma, nuval, 'o', label=r'$\nu(T)$',        linewidth=1.5)

if includefiles:
	raw = loadtxt('%s/PRE71_056408_Fig4.png.solid_line.dat' % subdir)
	ax.plot(raw[:,0],raw[:,1],'-', label=r'PRE71,056408, Fig.4, $T_{\rm i}=T_{\rm e}$', linewidth=1.0)
	raw = loadtxt('%s/PRE71_056408_Fig4.png.dotted_line.dat' % subdir)
	ax.plot(raw[:,0],raw[:,1],'--', label=r'PRE71,056408, Fig.4, $T_{\rm i}=1000{\rm K}$', linewidth=1.0)

	raw = loadtxt('%s/PRE71_056408_Fig4.png.filled_squares.dat' % subdir)
	ax.plot(raw[:,0],raw[:,1],linestyle='None', marker='s', label=r'PRE 71,056408,Fig.4, $T_{\rm i}=T_{\rm e}$, MD-sim', linewidth=1.0)
	raw = loadtxt('%s/PRE71_056408_Fig4.png.empty_squares.dat' % subdir)
	ax.plot(raw[:,0],raw[:,1],linestyle='None', marker='s', label=r'PRE 71,056408,Fig.4, $T_{\rm i}=1000{\rm K}$, MD-sim', linewidth=1.0, markerfacecolor='None')


ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='lower right')


###################################
	
plt.savefig('collision_frequency.pdf') # Must occure before show()
plt.show()
