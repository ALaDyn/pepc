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
from scipy import polyval
from numpy import mean, cov, corrcoef, var
from scipy.stats import t as student_tau

#print "Plots data from energy.dat and nose_hoover.dat"
#print "Usage %s [show] [useE0]" % sys.argv[0]

#locale.setlocale(locale.LC_ALL, "en_US.utf8") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
rc('font', family='serif', size=28)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

filename='energy.dat'


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
useE0    = 'useE0' in sys.argv

fig = plt.figure(figsize=(16,8))
#fig.suptitle(filename, fontsize=26)
fig.subplots_adjust(right=0.975, bottom=0.1,left=0.1,top=0.95,wspace=0.225,hspace=0.1) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
prop = matplotlib.font_manager.FontProperties(size=16)


####################################
ax=fig.add_subplot(1,2,1)
plt.xlabel(datanames[0])
plt.ylabel("Energy [Ry]")
plt.grid(True, which='both')
#ax.set_yscale('log')

for colidx in [1,4,5, 6,7,9]:
	plt.plot(data[:,0],data[:,colidx],label=datanames[colidx], linewidth=1.5,color=cm.hsv(24*colidx))#, markersize=1.0, markeredgewidth=1.0,marker=markers[colidx])

plt.legend(numpoints=1, prop=prop, ncol=2, frameon=True, loc='upper right')#, bbox_to_anchor=(1.05, 1), borderaxespad=0.)


####################################

lims=plt.axis()
#plt.gca().set_ylim([-0.6,3.4 ])
#plt.gca().set_xlim([ 0. ,5.75])
#plt.xticks(arange(lims[0],lims[1],2.0)) 
plt.minorticks_on()


####################################
####################################

filename='nose_hoover.dat'
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

abnm               = 5.2917720859E-2
kB                 = 1.
eps0               = 1./(4.*pi)
unit_t0_in_s       =       4.8377687E-17
unit_t0_in_fs      = unit_t0_in_s / 1.E-15
unit_qe            = -sqrt(2.)
eV_Ry              = 13.60569193
unit_me            = 1./2.


if 'ne_nm3' in simparams.keys():
  simparams['ne'] = simparams['ne_nm3'] * abnm**3
  
simparams['wpl'] = sqrt(16*pi*simparams['ne'] )

# we skip everything before switching on the laser
for startidx in range(len(data[:,0])):
  if data[startidx,0] > simparams['t0']:
    break
    
    
###################################
# linear fit for He+Hi  wit   y = a*x + b
FITCOL = 15
fitrawx=data[startidx:-1,0]
fitrawy=data[startidx:-1, FITCOL]

def moving_avg_y(xdata, ydata, window):
  weights = repeat(1.0, window) / window
  return xdata[window-1:-window+1], convolve(ydata, weights, mode='same')[window-1:-window+1]

delt    = fitrawx[1]-fitrawx[0]
Tlaser  = 2*pi/(simparams['w_over_wpl']*simparams['wpl'])*unit_t0_in_fs
nTLaser = Tlaser/delt

# moving average over 2..4 laser periods to reduce artificial statistical fluctuations
fitrawx, fitrawy = moving_avg_y(fitrawx, fitrawy, 4*nTLaser)

mycov = cov(fitrawx, fitrawy,ddof=1)
a     = mycov[0,1]/mycov[0,0]
b     = mean(fitrawy) - a*mean(fitrawx)
yfit  = polyval( [a,b], data[startidx:-1,0])

nobserv = len(fitrawx)
# taken from http://mathworld.wolfram.com/LeastSquaresFitting.html
s = sqrt((mycov[1,1] - a*mycov[0,1])/(nobserv -2))
a_std   = s*sqrt(1./nobserv + mean(fitrawx)**2/mycov[0,0])
b_std   = s/sqrt(mycov[0,0])
CONFIDENCE_LEVEL=0.95

a_interval = student_tau.interval(CONFIDENCE_LEVEL, nobserv, loc=a, scale=a_std)
b_interval = student_tau.interval(CONFIDENCE_LEVEL, nobserv, loc=b, scale=b_std)
c_interval = student_tau.interval(CONFIDENCE_LEVEL, nobserv, loc=0, scale=1.0)

#print "Fit result:           %f * t[fs] + %f" % (a, b)
#print "confidence_intervals: ", c_interval
#print "a_interval:           ", a_interval
#print "b_interval:           ", b_interval

simparams['fit_slope']            = a
simparams['fit_interval_slope']   = a_interval
simparams['fit_interval_factors'] = c_interval
simparams['fit_rel_error_slope']  = (c_interval[1] * a_std)/a
simparams['Te_Ryd'] = simparams['Te'] / eV_Ry
simparams['Ti_Ryd'] = simparams['Ti'] / eV_Ry

# function for temp-dependent collision frequency
def nu(Temp, dU_dT, s):
        if ('v0_over_vtherm' in s.keys()) and not useE0:
          # Use expression for fixed v0/vtherm
          Etherm = 3./2.*kB*Temp
          return 1./s['v0_over_vtherm']**2 / Etherm * dU_dT
        else:
          # Use expression for fixed E0
          return s['w_over_wpl']**2 * 2./(eps0*s['E0']**2) * s['ne'] * dU_dT

simparams

simparams['Gamma'] =     simparams['Z']*unit_qe**2/(kB*simparams['Te_Ryd'])*(4.*pi/3.*simparams['ne'])**(1./3.)
simparams['Theta'] = (kB*simparams['Te_Ryd'])*(3*pi*pi*simparams['ne'])**(-2./3.)

simparams['nuval']     = nu(simparams['Te_Ryd'], a*unit_t0_in_fs/(simparams["Ne"] + simparams["Ni"]), simparams)
simparams['nuval_max'] = simparams['nuval'] * (1.+simparams['fit_rel_error_slope'])
simparams['nuval_min'] = simparams['nuval'] * (1.-simparams['fit_rel_error_slope'])

simparams['nuval_wpl']     = simparams['nuval']      / simparams['wpl']
simparams['nuval_wpl_max'] = simparams['nuval_max']  / simparams['wpl']
simparams['nuval_wpl_min'] = simparams['nuval_min']  / simparams['wpl']

simparams['w']     = simparams['w_over_wpl'] * simparams['wpl']
simparams['a_ii']  = (4.*pi/3. * simparams['ne'])**(-1./3.)
simparams['vth']   = sqrt(3.*kB*simparams['Te_Ryd']/unit_me)
simparams['vosc']  = simparams['v0_over_vtherm']*simparams['vth']
simparams['xosc']  = simparams['vosc'] / simparams['w']
simparams['xosc_over_a_ii']  = simparams['xosc'] / simparams['a_ii']


print "***************** NU/WPL: [%g ... %g ... %g]" % (simparams['nuval_wpl_min'], simparams['nuval_wpl'], simparams['nuval_wpl_max'])

# Dump parameters and results
f = open('./simresults.dat', 'w')
json.dump(simparams, f)
f.close()


###################################
ax = fig.add_subplot(1,2,2)
ax.set_xlabel(datanames[0])
ax.set_ylabel(r'$\frac{1}{N}\mathcal{H}_N\,\mathrm{[Ry]}$')
ax.grid(True, which='both')
#ax.set_yscale('log')
#ax.set_xlim([0.,30.])

for p in [7, 14, 15]:
  ax.plot(data[startidx:-1,0], data[startidx:-1,p]/(simparams["Ne"] + simparams["Ni"]),  '-', label=datanames[p], linewidth=1.5)

ax.plot(fitrawx, fitrawy/(simparams["Ne"] + simparams["Ni"]),  '-', label=datanames[15] + ' moving avg', linewidth=1.5)


ax.plot(data[startidx:-1,0], yfit/(simparams["Ne"] + simparams["Ni"]),  '--', label="%s - fit" % datanames[FITCOL], linewidth=0.5, color='black')

ax.legend(numpoints=1, prop=prop, ncol=2, frameon=True, loc='upper left')


####################################
####################################
	
plt.savefig('energy_nose_hoover.pdf') # Must occure before show()
if showplot:
  plt.show()
