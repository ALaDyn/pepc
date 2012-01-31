#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

numruns = 37

avg        = []
avgP       = []
avgCluster = []
avgEnergy  = []

avg        = 0
avgP       = 0
avgCluster = 0
avgEnergy  = 0

for i in range(1,numruns+1):
  print "Processing step %06d" % i

  filename        = "rundir-%06d/momentum_electrons_Kt.dat" % i
  filenameP       = "rundir-%06d/momentum_electrons.dat" % i
  filenameCluster = "rundir-%06d/cluster.dat" % i
  filenameEnergy  = "rundir-%06d/energy.dat" % i
 
  rawdata         = np.loadtxt(filename)
  rawdataP        = np.loadtxt(filenameP)
  rawdataCluster  = np.loadtxt(filenameCluster)
  rawdataEnergy   = np.loadtxt(filenameEnergy)
  
  avg        = avg        + rawdata
  avgP       = avgP       + rawdataP
  avgCluster = avgCluster + rawdataCluster
  avgEnergy  = avgEnergy  + rawdataEnergy
  
avg        = avg        / numruns
avgP       = avgP       / numruns
avgCluster = avgCluster / numruns
avgEnergy  = avgEnergy  / numruns


np.savetxt("./momentum_electrons_Kt_avg.dat", avg)
np.savetxt("./momentum_electrons_avg.dat",    avgP)
np.savetxt("./cluster_avg.dat",               avgCluster)
np.savetxt("./energy_avg.dat",                avgEnergy)

# ======== plotting part ============

plt.plot(avg[:,0], avg[:,1])
plt.xlabel("t / fs")
plt.ylabel("K(t)")

plt.savefig("./momentum_electrons_Kt_avg.pdf")

plt.show()
