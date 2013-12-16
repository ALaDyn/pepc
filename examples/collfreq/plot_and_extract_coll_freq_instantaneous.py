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
from scipy.interpolate import UnivariateSpline

#print "Plots data from energy.dat for instantaneous heating rate measurement"
#print "Usage %s [show]" % sys.argv[0]

#locale.setlocale(locale.LC_ALL, "en_US.utf8") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
rc('font', family='serif', size=28)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

filename='energy.dat'
fileout='simresults_data.dat'


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
          
showplot = 'show' in sys.argv

fig = plt.figure(figsize=(16,8))
#fig.suptitle(filename, fontsize=20)
fig.subplots_adjust(right=0.95, top=0.9, bottom=0.15,left=0.1, wspace=0.25) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
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


##################################
abnm = 5.2917720859E-2
kB   = 1.
eps0 = 1./(4.*pi)
unit_t0_in_s       =       4.8377687E-17
unit_t0_in_fs      = unit_t0_in_s / 1.E-15
unit_qe            = -sqrt(2.)
eV_Ry              = 13.60569193
unit_me            = 1./2.
unit_hbar          = 1.

if 'ne_nm3' in simparams.keys():
  simparams['ne'] = simparams['ne_nm3'] * abnm**3

simparams['wpl'] = sqrt(16*pi*simparams['ne'] )

simparams['Te_Ryd'] = simparams['Te'] / eV_Ry
simparams['Ti_Ryd'] = simparams['Ti'] / eV_Ry

simparams['Gamma'] =     simparams['Z']*unit_qe**2/(kB*simparams['Te_Ryd'])*(4.*pi/3.*simparams['ne'])**(1./3.)
simparams['Theta'] = (kB*simparams['Te_Ryd'])/(unit_hbar**2*(3*pi*pi*simparams['ne'])**(2./3.) / (2*unit_me))

simparams['data'] = fileout

# Dump parameters and results
f = open('./simresults.dat', 'w')
json.dump(simparams, f)
f.close()

    
#### perform fit for Etherm_electrons
# Fitfunction for thermal energy, time in simulation units
def Efunc(p, t):
	return p(t)
# derivative of fitfunction, time in simulation units
def dEfunc(p, t):
        return p(t,nu=1)
# function for time-dependent collision frequency, time in simulation units
def nu(p, t, s):
        if ('v0_over_vtherm' in s.keys()):
          # Use expression for fixed v0/vtherm
          Etherm = Efunc(p, t)
          return 1./s['v0_over_vtherm']**2 / Etherm * dEfunc(p, t)
        else:
          # Use expression for fixed E0
          return s['w_over_wpl']**2 * 2./(eps0*s['E0']**2) * s['ne'] * dEfunc(p, t)

Tlaser  = 2*pi/(simparams['w_over_wpl']*simparams['wpl'])*unit_t0_in_fs
SKIPCYCLES=5

# we skip everything before switching on the laser and the first SKIPCYCLE cycles
for startidx in range(len(data[:,0])):
  if data[startidx,0] > simparams['t0'] + max(SKIPCYCLES*Tlaser, 3.):
    break
    
polydegree = 3
smoothing  = 1

#polydegree = 2
#smoothing  = 12.5*max(simparams['Te'], 0.0) # probably a good rule of thumb - who knows :-)


xvals     = data[startidx:-1,0] / unit_t0_in_fs #scaled to simulation time units - otherwise, the derivative does not work as intended
yvals     = data[startidx:-1,7]

p1 = UnivariateSpline(xvals, yvals, k=polydegree, s=smoothing)

yvals_fit  =  Efunc(p1, xvals)
dyvals_fit = dEfunc(p1, xvals)

for startidx in range(len(data[:,0])):
  if dyvals_fit[startidx] > 1e-4:
    break

xvals     =  xvals[startidx:-1]
yvals_fit =  yvals_fit[startidx:-1]
dyvals_fit= dyvals_fit[startidx:-1]

nuvals     = nu(   p1, xvals, simparams)
nuvals_wpl = nuvals  / simparams['wpl']
Tvals      = Efunc(p1, xvals) * 2./(3.*kB) * eV_Ry
Gammavals  = simparams['Z']*unit_qe**2/(kB*Tvals/eV_Ry)*(4.*pi/3.*simparams['ne'])**(1./3.)
Thetavals  = (kB*Tvals/eV_Ry)*(3*pi*pi*simparams['ne'])**(-2./3.)
nevals     = ones(Tvals.shape) * simparams['ne']

rawdataout = array([nevals, Tvals, Gammavals, Thetavals, nuvals, nuvals_wpl]).transpose()
savetxt(fileout, rawdataout)


###################################
ax = fig.add_subplot(1,2,1)
ax.set_xlabel(datanames[0])
ax.set_ylabel("Energy [Ry]")
ax.grid(True, which='both')
#ax.set_yscale('log')
#ax.set_xlim([0.,30.])

ax.plot(data[:,0], data[:,7],  '-', label=datanames[7], linewidth=1.5)
ax.plot(xvals*unit_t0_in_fs, yvals_fit,            '-', label='spline fit', linewidth=1.5)

ax.legend(numpoints=1, prop=prop, ncol=1, frameon=True, loc='lower right')
plt.minorticks_on()

####################################
ax = fig.add_subplot(1,2,2)
ax.set_xlabel(r'$\Gamma$')
ax.set_ylabel(r'$\nu_{\rm ei} / \omega_{\rm pl}$')
ax.grid(True, which='both')
minG  = min(Gammavals)
maxG  = max(Gammavals)
#ax.set_xlim(array([minG, maxG]))
ax.set_xscale('log')
ax.set_yscale('log')


ax.plot(Gammavals, nuvals_wpl,     '-', label=r'$\nu(T)$',        linewidth=1.5)

####################################
plt.savefig('collision_frequency_instantaneous.pdf') # Must occure before show()
if showplot:
  plt.show()
