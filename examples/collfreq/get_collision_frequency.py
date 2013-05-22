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
from scipy import optimize
from scipy import odr

locale.setlocale(locale.LC_ALL, "en_US") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
rc('font',   family='serif', size=24)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

filename='energy.dat'

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
           '$U_{\mathrm{pot}}^{\mathrm{(total)}}$',
           '$U_{\mathrm{pot}}^{\mathrm{(near field)}}$',
           '$U_{\mathrm{pot}}^{\mathrm{(far field)}}$',
           '$U_{\mathrm{kin}}^{\mathrm{(electrons)}}$',
           '$U_{\mathrm{kin}}^{\mathrm{(ions)}}$',
           '$U_{\mathrm{kin}}^{\mathrm{(total)}}$',
           '$U_{\mathrm{kin,w/o\ drift}}^{\mathrm{(electrons)}}$',
           '$U_{\mathrm{kin,w/o\ drift}}^{\mathrm{(ions)}}$',
           '$U^{\mathrm{(tot)}}$',
           '$T_{\mathrm{e}}^0$',
           '$T_{\mathrm{e}}^{\mathrm{uncor}}$',
           '$\chi_{\mathrm{e}}$',
           '$\Delta T_{\mathrm{e}}$',
           '$T_{\mathrm{i}}^0$',
           '$T_{\mathrm{i}}^{\mathrm{uncor}}$',
           '$\chi_{\mathrm{i}}$',
           '$\Delta T_{\mathrm{i}}$'
          ]

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
  simparams = {'defaultfile' : True, 't0': 3.5, 'w_over_wpl' : 3.0, 'E0' : 1.56542, 'wpl' : 0.272921, 'Z' : 1, 'ne' : 0.00148185, 'v0' : 5.40775}
  f = open('./simparams.dat', 'w')
  json.dump(simparams, f)
  f.close()

if simparams['defaultfile']:
  print ">>>>>>>>>> simparams.dat contains default vaules - you should correct this and set defaultfile=False therein before continuing <<<<<<<<<<"
  sys.exit(0)


#### perform fit for Etherm_electrons
kB   = 1.
eps0 = 1./(4.*pi)
unit_t0_in_s       =       4.8377687E-17
unit_t0_in_fs      = unit_t0_in_s / 1.E-15
unit_qe            = -sqrt(2.) 

# Fitfunction for thermal energy, time in simulation units
def Efunc(p, t):
	return p[0] * (t-p[2])**p[1] + p[3]
        
# derivative of fitfunction, time in simulation units
def dEfunc(p, t):
	return p[0]*p[1] * (t-p[2])**(p[1]-1)
        
# function for time-dependent collision frequency, time in simulation units
def nu(p, t, s):
        if 'v0_over_vtherm' in s.keys():
          # Use expression for fixed v0/vtherm
          Etherm = Efunc(p, t)
          return 1./s['v0_over_vtherm'] * 4. / Etherm * dEfunc(p, t) / s['wpl']
        else:
          # Use expression for fixed E0
          return s['w_over_wpl']**2 * 2./(eps0*s['E0']**2) * s['ne'] * dEfunc(p, t) / s['wpl']

# we skip everything before switching on the laser
for startidx in range(len(data[:,0])):
  if data[startidx,0] > simparams['t0']:
    break

xvals = data[startidx:-1,0] / unit_t0_in_fs #scaled to simulation time units - otherwise, the derivative does not work as intended
yvals = data[startidx:-1,7]

p0 = [1., 1., 0., 0.]
# actually perform the fitting process
o = odr.ODR(odr.RealData(xvals, yvals), odr.Model(Efunc), p0, maxit=10000)
o.set_job(fit_type=2)
#o.set_iprint(init=2, iter=2, final=2)
out = o.run()
p1  = out.beta

yvals_fit = Efunc(p1, xvals)

out.pprint()


###################################
ax = fig.add_subplot(2,2,1)
ax.set_xlabel(datanames[0])
ax.set_ylabel("Energy [Ry]")
ax.grid(True, which='both')
#ax.set_yscale('log')
ax.set_xlim([0.,30.])

ax.plot(xvals*unit_t0_in_fs, yvals,    '-', label=datanames[7],             linewidth=1.5)
ax.plot(xvals*unit_t0_in_fs, yvals_fit,'-', label='fit: %s' % datanames[7], linewidth=1.5)

if includefiles:
	raw = loadtxt('%s/PRE71_056408_Fig3.png.dat' % subdir)
	ax.plot(raw[:,0],raw[:,1],'--', label='PRE71,056408,Fig.3', linewidth=1.0)

ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='lower right')
plt.minorticks_on()

###################################
ax = fig.add_subplot(2,2,3)
ax.set_xlabel(datanames[0])
ax.set_ylabel(r'$\nu_{\rm ei} / \omega_{\rm pl}$')
ax.grid(True, which='both')
ax.set_yscale('log')
ax.set_xlim([0.,30.])

nuvals     = nu(    p1, xvals, simparams)

ax.plot(xvals*unit_t0_in_fs,nuvals,     '-', label=r'$\nu(t)$',        linewidth=1.5)

ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='upper right')

####################################
ax = fig.add_subplot(2,2,2)
ax.set_xlabel(r'T [Ry]')
ax.set_ylabel(r'$\nu_{\rm ei} / \omega_{\rm pl}$')
ax.grid(True, which='both')
minT = min(yvals_fit) * 2./(3.*kB)
maxT = max(yvals_fit) * 2./(3.*kB)
#ax.set_xlim(array([minT, maxT]))
#ax.set_xscale('log')
ax.set_yscale('log')

Tvals = Efunc(p1, xvals) * 2./(3.*kB)
ax.plot(Tvals,nuvals,     '-', label=r'$\nu(T)$',        linewidth=1.5)

ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='lower left')

####################################
ax = fig.add_subplot(2,2,4)
ax.set_xlabel(r'$\Gamma$')
ax.set_ylabel(r'$\nu_{\rm ei} / \omega_{\rm pl}$')
ax.grid(True, which='both')
Gamma = simparams['Z']*unit_qe**2/(kB*Tvals)*(4.*pi/3.*simparams['ne'])**(1./3.)
minG  = min(Gamma)
maxG  = max(Gamma)
#ax.set_xlim(array([minG, maxG]))
ax.set_xscale('log')
ax.set_yscale('log')


ax.plot(Gamma, nuvals,     '-', label=r'$\nu(T)$',        linewidth=1.5)

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

####################################

	
plt.savefig('collision_frequency.pdf') # Must occure before show()
plt.show()
