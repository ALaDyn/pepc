#!/usr/bin/python

import matplotlib.pyplot as plt
from matplotlib import rc
from numpy import *
import matplotlib.font_manager
import matplotlib.cm as cm
import matplotlib as mpl
import locale
import os
import sys
import json

#print "Plots data from simresults.dat with nose_hoover thermostat"
#print "Usage: %s [show]" % sys.argv[0]

DATADIR='../data_from_others/'

def mycm(k):
  bla = ['k', 'r', 'b', 'g', 'm', 'orange'] # mpl.rcParams['axes.color_cycle']
  return bla[k % len(bla)]

#locale.setlocale(locale.LC_ALL, "en_US.utf8") # fuer tausendertrennzeichen im format

#rc('text', usetex=True)
rc('font', family='serif', size=24)
rc('legend', fontsize='small')
rc('xtick', labelsize='small')
rc('ytick', labelsize='small')

rc('xtick.major', size='8')
rc('xtick.minor', size='4')
rc('ytick.major', size='8')
rc('ytick.minor', size='4')


legendfontsize=14

filename='simresults.dat'
omitfile="_omit_directory"

rawdata=[]
rawdata_int=[]

def get_immediate_subdirectories(dir):
    return [name for name in os.listdir(dir)
            if os.path.isdir(os.path.join(dir, name))]

for dir in get_immediate_subdirectories('.'):

  currfilename = "./%s/%s" % (dir, filename)
  omitfilename = "./%s/%s" % (dir, omitfile)

  if os.path.exists(currfilename) and not os.path.exists(omitfilename):
    f = open(currfilename, 'r')
    if dir[0:2] == 'Te':
      rawdata.append(json.load(f))
      rawdata[-1]['dir'] = dir
    elif 'int' in dir:
      rawdata_int.append(json.load(f))
      rawdata_int[-1]['dir'] = dir

    f.close()
    
print "Loaded %d datasets" % len(rawdata) 

###############################################################################

columnnames = {'Ni'            : r'$N_\mathrm{ion}$',
               'Ne'            : r'$N_\mathrm{el }$',
               'nuval'         : r'$\nu_\mathrm{ei}$',
               'nuval_min'     : r'$\nu_\mathrm{ei}(\mathrm{min})$', 
               'nuval_max'     : r'$\nu_\mathrm{ei}(\mathrm{max})$',
               'nuval_wpl'     : r'$\nu_\mathrm{ei} / \omega_\mathrm{pl}$',
               'nuval_wpl_min' : r'$\nu_\mathrm{ei}(\mathrm{min}) / \omega_\mathrm{pl}$', 
               'nuval_wpl_max' : r'$\nu_\mathrm{ei}(\mathrm{max}) / \omega_\mathrm{pl}$',
               'Ti'            : r'$T_\mathrm{ion}\,\mathrm{[eV]}$',
               'Te'            : r'$T_\mathrm{el }\,\mathrm{[eV]}$',
               'Theta'         : r'$\Theta$',
               'Gamma'         : r'$\Gamma$',
               'ne'            : r'n_\mathrm{el}\,\mathrm{[a_\mathrm{B}^{-3}]}$', 
               'w_over_wpl'    : r'$\omega / \omega_\mathrm{pl}$',
               'v0_over_vtherm': r'$v_\mathrm{osc} / v_\mathrm{therm}$',
               'wpl'           : r'$\omega_\mathrm{pl}$', 
               'ne_nm3'        : r'$n_\mathrm{el}\,\mathrm{[nm^{-3}]}$', 
               'Z'             : r'$Z$',
               't0'            : r'$t_0$',
               'Ti_Ryd'        : r'$T_\mathrm{ion}\,\mathrm{[Ryd]}$',
               'Te_Ryd'        : r'$T_\mathrm{el }\,\mathrm{[Ryd]}$',
               'xosc_over_a_ii': r'$x_\mathrm{osc} / d_\mathrm{ii}$',
               }
markers = ['*','<','>','D','H','^','d','h','o','p','s','v','x','|','None',' ','']

columnint = { 'ne'         : 0,
              'Te'         : 1,
              'Gamma'      : 2,
              'Theta'      : 3,
              'nuval'      : 4,
              'nuval_wpl'  : 5
            }

###############################################################################

plotcounter = 1
MINY=1.e-6

def buildnewfig(xlab, ylab):
  global plotcounter
  # prepare plot
  fig = plt.figure(figsize=(12,8))
  fig.subplots_adjust(right=0.975, bottom=0.125,left=0.125,top=0.95,wspace=0.225,hspace=0.1) # http://matplotlib.sourceforge.net/faq/howto_faq.html#move-the-edge-of-an-axes-to-make-room-for-tick-labels
  ax = fig.add_subplot(1,1,1)
  plt.xlabel(xlab)
  plt.ylabel(ylab)
  #plt.grid(True, which='both')
  plotcounter = 0
  ax.set_xscale('log')
  ax.set_yscale('log')

def do_plot(xcol, ycol, ycol_max, ycol_min, label, inclusions, newfig, showerror, **kwargs):
  global plotcounter
  
  # extract data
  xvals     = []
  yvals     = []
  yerr_max  = []
  yerr_min  = []
  
  print "Plotting %s(%s) with filter" % (ycol, xcol), inclusions
  
  for line in rawdata:
    includeline = True
    
    for k in inclusions.keys():
      if line[k] != inclusions[k]:
        includeline = False
    
    if (includeline):
      xvals.append(line[xcol])
      yvals.append(line[ycol])
      
      if showerror:
        yerr_max.append(line[ycol_max] - line[ycol])
        ylower = maximum(MINY, line[ycol_min]) # see https://github.com/matplotlib/matplotlib/issues/163
        yerr_min.append(line[ycol] - ylower)
  
  if len(xvals) < 1:
    print "No data found - your filter might be invalid"
    
  # sort data with respect to x values
  xvals = array(xvals)
  yvals = array(yvals)
  yerr_max = array(yerr_max)
  yerr_min = array(yerr_min)
  sortidx = argsort(xvals)
  #print sortidx
  xvals = xvals[sortidx]
  yvals = yvals[sortidx]
  if showerror:
    yerr_max = yerr_max[sortidx]
    yerr_min = yerr_min[sortidx]

  if newfig:
    buildnewfig(columnnames[xcol], columnnames[ycol])
  
  ax = plt.gca()
  
  if showerror:
    ax.errorbar(xvals, yvals, yerr=[yerr_min, yerr_max], label=label, 
            fmt='%s' % markers[plotcounter], markersize = 10.0, capsize=10.0, elinewidth=2.0, color=mycm(plotcounter), linewidth=1.)
  else:
    ax.plot(xvals, yvals, label=label, linewidth=0.5, 
            marker=markers[plotcounter], markersize = 10.0, color=mycm(plotcounter), **kwargs )
            
  plotcounter = plotcounter + 1


def do_plot_int(xcol, ycol, label, inclusions, newfig):
  global plotcounter
  
  founddata = False
  print "Plotting instantaneous data for %s(%s) with filter" % (ycol, xcol), inclusions
  
  if newfig:
    buildnewfig(columnnames[xcol], columnnames[ycol])

  ax = plt.gca()
  dolabel = True
  
  for line in rawdata_int:
    includeline = True
    
    for k in inclusions.keys():
      if line[k] != inclusions[k]:
        includeline = False
    
    if (includeline):
      founddata = True
      datafilename = "%s/%s" % (line['dir'], line['data'])
      realraw = loadtxt(datafilename)
      
      if (dolabel):
        dolabel = False
        ax.plot(realraw[:,columnint[xcol]], realraw[:,columnint[ycol]], color=mycm(plotcounter), label=label)#, linewidth=1.)
      else:
        ax.plot(realraw[:,columnint[xcol]], realraw[:,columnint[ycol]], color=mycm(plotcounter))#, label=("%s = %f" % (columnnames[label], labelval)), linewidth=1.)
            
  if founddata:
    plotcounter = plotcounter + 1


def do_data_plot(datafile, newfig, showerror, **kwargs):
  global plotcounter
  
  # load data
  filename = "%s/%s" % (DATADIR, datafile)
  
  rawdata = loadtxt(filename)
  
  if newfig:
    buildnewfig('', '')
  
  ax = plt.gca()
  
  if 'markersize' in kwargs.keys():
    kwargs['marker'] = markers[plotcounter]
    
  if 'datafactor' in kwargs.keys():
    factor =  kwargs.pop('datafactor')
  else:
    factor = 1.0
    
  if 'showlabel' in kwargs.keys():
    if kwargs.pop('showlabel'):
     kwargs.pop('label')
  
  if showerror:
    ax.errorbar(xvals, factor*yvals, yerr=[yerr_min, yerr_max], 
            fmt=markers[plotcounter], markersize = 7.5, capsize=7.5, elinewidth=1.5, **kwargs) #, color=mycm(plotcounter)
  else:
    ax.plot(rawdata[:,0], factor*rawdata[:,1], color=mycm(plotcounter), markerfacecolor='None', **kwargs) # 
            
  plotcounter = plotcounter + 1



###################################
        
do_data_plot('PhysRevLett.103.065005_Fig1_Hilse.dat',
             True,
             False,
             label=r'Hilse et al.',
             linewidth=0., markersize = 10
             )        

do_data_plot('PhysRevLett.103.065005_Fig1_PfalznerGibbon.dat',
             False,
             False,
             label=r'Pfalzner & Gibbon',
             linewidth=0., markersize = 10
             )        

do_data_plot('PhysRevLett.103.065005_Fig1_GrinenkoGericke.dat',
             False,
             False,
             label=r'Grinenko & Gericke',
             linewidth=2, linestyle='--'
             )        

do_data_plot('PhysRevLett.103.065005_Fig1_GrinenkoGerickePolarization.dat',
             False,
             False,
             label=r'Grinenko & Gericke, polarization only',
             linewidth=2, linestyle='-.'
             )        

do_data_plot('PhysRevLett.103.065005_Fig1_Decker.dat',
             False,
             False,
             label=r'Decker et al.',
             linewidth=2, linestyle='--'
             )        

do_data_plot('PhysRevLett.103.065005_Fig1_CaubleRozmus.dat',
             False,
             False,
             label=r'Cauble & Rozmus',
             linewidth=0., markersize = 10
             )        

do_data_plot('PhysRevLett.103.065005_Fig1_Bornath.dat',
             False,
             False,
             label=r'Bornath et al.',
             linewidth=2, linestyle='--'
             )        

do_plot('Gamma', 
        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
        'this work', 
        {'v0_over_vtherm' :  0.2,
         'w_over_wpl'     :  3.0,
         'ne_nm3'         : 10.0},
        False,
        True)

plt.xlabel(columnnames['Gamma'])
plt.ylabel(columnnames['nuval_wpl'])

plt.legend(loc='upper left', ncol=1,prop={'size':legendfontsize})
plt.gca().set_xlim([1e-2, 2e1])
plt.gca().set_ylim([2e-3, 1e1])


plt.savefig('nu_over_wpl(Gamma)_with_data.pdf') # Must occure before show()

###################################

do_plot_int('Gamma', 
        'nuval_wpl',
        'instantaneous measurement', 
        {'v0_over_vtherm' :  0.2,
         'w_over_wpl'     :  3.0,
         'ne_nm3'         : 10.0},
        True)

do_plot('Gamma', 
        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
        'Nose-Hoover thermostat', 
        {'v0_over_vtherm' :  0.2,
         'w_over_wpl'     :  3.0,
         'ne_nm3'         : 10.0},
        False,
        True)

plt.legend(loc='upper left', ncol=1,prop={'size':legendfontsize})
plt.gca().set_xlim([1e-2, 2e1])
plt.gca().set_ylim([1e-2, 1e0])


plt.savefig('nu_over_wpl(Gamma)_with_instantaneous.pdf') # Must occure before show()

###################################

#newplot = True
#for dens in [0.001, 0.1, 10.0]:
#  do_plot('Gamma', 
#        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
#        '%s = %g' % (columnnames['ne_nm3'], dens),
#        {'v0_over_vtherm' :  0.2,
#         'w_over_wpl'     :  3.0,
#         'ne_nm3'         : dens},
#        newplot,
#	True)
#
#  newplot=False
#
#plt.legend(loc='lower right', ncol=1, prop={'size':legendfontsize})
#plt.gca().set_xlim([1e-3, 1e2])
#plt.gca().set_ylim([2e-3, 1e1])
#
#plt.savefig('nu_over_wpl(Gamma,ne)_with_data.pdf') # Must occure before show()

###################################

#do_plot_int('Theta', 
#        'nuval_wpl',
#        '%s = %g, inst. meas' % (columnnames['ne_nm3'], 10.0),
#         {'v0_over_vtherm' :  0.2,
#         'w_over_wpl'     :  3.0,
#         'ne_nm3'         : 10.0},
#        True)
#
#for dens in [0.001, 0.1, 10.0]:
#  do_plot('Theta', 
#        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
#        '%s = %g' % (columnnames['ne_nm3'], dens),
#        {'v0_over_vtherm' :  0.2,
#         'w_over_wpl'     :  3.0,
#         'ne_nm3'         : dens},
#        False,
#	True)
#
#
#plt.legend(loc='lower left', ncol=1, prop={'size':legendfontsize})
#plt.savefig('nu_over_wpl(Theta)_with_data.pdf') # Must occure before show()

###################################

mytvals=[['2'  , 2.0, 0.00725/2e-5],
         ['7.7', 7.7, 0.00055/1.06e-6],
         #['33' ,33.0],
         ]
         
for tval in mytvals:
  newplot=True

  for vv in [0.2]:#, 2.0]:
    do_plot('w_over_wpl', 
        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
        r'%s = %g, $v_\mathrm{osc} / v_\mathrm{therm} = %g$' % (columnnames['Te'], tval[1], vv),
        {'v0_over_vtherm' :  vv,
         'Te'             :  tval[1],
         'ne_nm3'         : 10.0},
        newplot,
        True)
    newplot=False
  
  myx = array([1e-4, 1e4])
  myy =tval[2]* myx**(-7./2.)
  
  plt.gca().plot(myx, myy, color='gray', linewidth=1.0)
  
  if tval[1] != 7.7:
  
    do_data_plot('../data_rostock/Rostock_Te%seV_ne1.0E22percm3_Re(nu_Born).dat' % tval[0],
             False,
             False,
             label='1st Born',
             linewidth=1.5, linestyle='--'
             )
             
    plotcounter = plotcounter - 1
    do_data_plot('../data_rostock/Rostock_Te%seV_ne1.0E22percm3_Re(nu_LB).dat' % tval[0],
             False,
             False,
             label='LB (dyn.scr.Born)',
             linewidth=1.5, linestyle='-.'
             )        

    plotcounter = plotcounter - 1
    do_data_plot('../data_rostock/Rostock_Te%seV_ne1.0E22percm3_Re(nu_TM).dat' % tval[0],
             False,
             False,
             label='T-Matrix(ladder approx.)',
             linewidth=1.5, linestyle=':'
             )        

    plotcounter = plotcounter - 1
    do_data_plot('../data_rostock/Rostock_Te%seV_ne1.0E22percm3_Re(nu_GDW).dat' % tval[0],
             False,
             False,
             label='Gould-deWitt',
             linewidth=2.25, linestyle='-'
             )      
  else:  
    do_data_plot('../data_rostock/Rostock_Te%seV_ne1.0E22percm3_Re(nu_TM).dat' % tval[0],
             False,
             False,
             label='Gould-deWitt',
             linewidth=2.25, linestyle='-'
             )        


  plt.gca().set_xlim([0.05, 200])
  plt.gca().set_ylim([1e-4, 1.0])


  plt.legend(loc='lower left', ncol=1, prop={'size':legendfontsize})
  plt.savefig('nu_over_wpl(w_over_wpl,T%seV)_with_data.pdf' % tval[0]) # Must occure before show()

###################################

mytvals=[ [ 2.0,  '2.0', 10.0, '10.0', 45.67/3.769],
          [50.0, '50.0', 10.0, '10.0', 11.67/7.002],
	      #[50.0, '50.0',  0.1 , '0.1', 1.0],
           ]

for do_rescale in [False, True]:
 newplot=True
 for tval in mytvals:
  do_plot('v0_over_vtherm', 
        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
        '%s=%g, %s=%g' % (columnnames['Te'], tval[0], columnnames['ne_nm3'], tval[2]),
        {'w_over_wpl'     :  3.0,
         'Te'             :   tval[0],
         'ne_nm3'         :   tval[2]},
        newplot,
        True)
  
  if do_rescale:
    datfac = tval[4]
  else:
    datfac = 1.
  
  plotcounter = plotcounter - 1
  do_data_plot('SovietPhysicsJETP20_1510_SilinAll_Te%s_ne%s.dat' % (tval[1], tval[3]),
             False,
             False,
             label=r'Silin',
             #label=r'Silin, eq. (3.10), $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
             linewidth=1.5, linestyle='-', datafactor=datfac, showlabel=newplot
             )        
  plotcounter = plotcounter - 1
  do_data_plot('SovietPhysicsJETP20_1510_SilinLoV_Te%s_ne%s.dat' % (tval[1], tval[3]),
             False,
             False,
             label=r'Dawson & Oberman',
             #label=r'Dawson & Oberman = Silin, eq. (3.12), $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
             linewidth=1.5, linestyle=':', datafactor=datfac, showlabel=newplot
             )        
  #plotcounter = plotcounter - 1
  #do_data_plot('SovietPhysicsJETP20_1510_SilinHiV_Te%s_ne%s.dat' % (tval[1], tval[3]),
  #           False,
  #           False,
  #           #label=r'Silin, eq. (3.13), $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
  #           linewidth=1.5, linestyle='--', datafactor=datfac, showlabel=newplot
  #           )        

  if do_rescale:
    plotcounter = plotcounter - 1
    do_data_plot('SovietPhysicsJETP20_1510_SilinLoVLangdon_Te%s_ne%s.dat' % (tval[1], tval[3]),
             False,
             False,
             label=r'Dawson & Oberman with Langdon',
             #label=r'Silin, eq. (3.13) with Langdon corr., $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
             linewidth=1.5, linestyle='-.', datafactor=datfac, showlabel=newplot
             )        
             
  if do_rescale:
    datfac = tval[4]/0.447
  else:
    datfac = 1.

  #plotcounter = plotcounter - 1
  #do_data_plot('SovietPhysicsJETP20_1510_SilinLangdon_Te%s_ne%s.dat' % (tval[1], tval[3]),
  #           False,
  #           False,
  #           label=r'Silin with Langdon',
  #           #label=r'Silin, eq. (3.13) with Langdon corr., $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
  #           linewidth=1.5, linestyle='--', datafactor=datfac, showlabel=newplot
  #           )        
             
  newplot = False
 plt.gca().set_ylim([1e-7, 100.0]) 
 plt.gca().set_xlim([0.05, 120.0]) 

 plt.legend(loc='upper center', ncol=2, prop={'size':legendfontsize})

 if do_rescale:
   plt.savefig('nu_over_wpl(v0_over_vtherm)_with_data.pdf') # Must occure before show()
 else:
   plt.savefig('nu_over_wpl(v0_over_vtherm)_with_data,no_rescale.pdf') # Must occure before show()

###################################

for do_rescale in [True]:
 newplot=True
 for tval in mytvals:
  do_plot('xosc_over_a_ii', 
        'nuval_wpl', 'nuval_wpl_max', 'nuval_wpl_min', 
        '%s=%g, %s=%g' % (columnnames['Te'], tval[0], columnnames['ne_nm3'], tval[2]),
        {'w_over_wpl'     :  3.0,
         'Te'             :   tval[0],
         'ne_nm3'         :   tval[2]},
        newplot,
        True)
  
  if do_rescale:
    datfac = tval[4]
  else:
    datfac = 1.

  plotcounter = plotcounter - 1
  do_data_plot('SovietPhysicsJETP20_1510_SilinAll_Te%s_ne%s_xosc.dat' % (tval[1], tval[3]),
             False,
             False,
             label=r'Silin',
             #label=r'Silin, eq. (3.10), $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
             linewidth=1.5, linestyle='-', datafactor=datfac, showlabel=newplot
             )        
  plotcounter = plotcounter - 1
  do_data_plot('SovietPhysicsJETP20_1510_SilinLoV_Te%s_ne%s_xosc.dat' % (tval[1], tval[3]),
             False,
             False,
             label=r'Dawson & Oberman',
             #label=r'Dawson & Oberman = Silin, eq. (3.12), $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
             linewidth=1.5, linestyle=':', datafactor=datfac, showlabel=newplot
             )        
  #plotcounter = plotcounter - 1
  #do_data_plot('SovietPhysicsJETP20_1510_SilinHiV_Te%s_ne%s_xosc.dat' % (tval[1], tval[3]),
  #           False,
  #           False,
  #           label=r'Silin w Langdon',
  #           #label=r'Silin, eq. (3.13) with Langdon corr., $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
  #           linewidth=1.5, linestyle='--', datafactor=datfac, showlabel=newplot
  #           )        

  if do_rescale:
    plotcounter = plotcounter - 1
    do_data_plot('SovietPhysicsJETP20_1510_SilinLoVLangdon_Te%s_ne%s_xosc.dat' % (tval[1], tval[3]),
             False,
             False,
             label=r'Dawson & Oberman with Langdon',
             #label=r'Silin, eq. (3.13) with Langdon corr., $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
             linewidth=1.5, linestyle='-.', datafactor=datfac, showlabel=newplot
             )        
             
  if do_rescale:
    datfac = tval[4]/0.447
  else:
    datfac = 1.

  #plotcounter = plotcounter - 1
  #do_data_plot('SovietPhysicsJETP20_1510_SilinLangdon_Te%s_ne%s_xosc.dat' % (tval[1], tval[3]),
  #           False,
  #           False,
  #           label=r'Silin w Langdon',
  #           #label=r'Silin, eq. (3.13) with Langdon corr., $T_\mathrm{e}=%g\,\mathrm{eV},\,n_\mathrm{e}=%g\,\mathrm{nm^{-3}}$' % (tval[0], tval[2]),
  #           linewidth=1.5, linestyle='--', datafactor=datfac, showlabel=newplot
  #           )        
             
  newplot = False

 plt.gca().set_ylim([1e-7,  50.0]) 
 plt.gca().set_xlim([1e-2, 120.0]) 
 
 plt.legend(loc='upper center', ncol=2, prop={'size':legendfontsize})
 if do_rescale:
   plt.savefig('nu_over_wpl(a0_over_dii)_with_data.pdf') # Must occure before show()
 else:
   plt.savefig('nu_over_wpl(a0_over_dii)_with_data,no_rescale.pdf') # Must occure before show()

###################################

#TODO: Plot fuer Gamma approx 10. und verschiedene Dichten / Temperaturen

if 'show' in sys.argv:
  plt.show()


