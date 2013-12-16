#!/lustre/jhome2/jzam04/jzam0415/programs/python/bin/python

from pylab import *
import json

clusterdata = np.loadtxt('./cluster_jackknife.dat')

simparams ={}

simparams['rion']     = clusterdata[-1,2]
simparams['rel']      = clusterdata[-1,3]
simparams['Nion']     = clusterdata[-1,4]
simparams['Nefree']   = clusterdata[-1,5]
simparams['Nebound']  = clusterdata[-1,6]

simparams['V_ionsphere']     = 4./3. * pi * simparams['rion']**3
simparams['ne_in_ionsphere'] = simparams['Nebound']/simparams['V_ionsphere']

simparams['V_esphere']       = 4./3. * pi * simparams['rel']**3
simparams['ne_in_esphere']   = simparams['Nebound']/simparams['V_esphere']

qe     = -sqrt(2.);
epsz   = 1./(4. * pi)
me     = 1./2.
t0infs = 4.8377687*10**(-17)/10**(-15)

simparams['wpl']  = sqrt(((simparams['ne_in_ionsphere'] * qe**2)/(epsz*me)))/t0infs
simparams['wMie'] = simparams['wpl'] / sqrt(3.)

simparams['wpl_ecloud']  = sqrt(((simparams['ne_in_esphere'] * qe**2)/(epsz*me)))/t0infs
simparams['wMie_ecloud'] = simparams['wpl_ecloud'] / sqrt(3.)

f = open('./simparams.dat', 'w')
json.dump(simparams, f)
f.close()
