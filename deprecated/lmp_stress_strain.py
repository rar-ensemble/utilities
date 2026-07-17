#!/usr/bin/python

import lammps_utils as lmp
import sys
import numpy as nmp
import os
import io_utils as io
import string
import math as m

if (len(sys.argv) < 2):
  print "Usage:  avg_lmp_windows [log.input]" 
  exit(0)

if (os.path.exists(sys.argv[1]) ):
  dat = lmp.read_log_dat(sys.argv[1])

  ntimes = nmp.shape(dat)[0]
  nvars = nmp.shape(dat)[1]

  inp = open(sys.argv[1],'r')
  tstep = 0
  thermo = 0
  while (inp):
    L = inp.readline().split()

    if (len(L) > 1 and L[0] == 'thermo'):
      thermo = int(L[1])

    if (len(L) > 1 and L[0] == 'timestep'):
      tstep = float(L[1])

    if (len(L) > 1 and L[0] == 'Step'):
      break

  aind = []
  pind = []
  for j in range(0,len(L)):
    
    if (L[j]=="Pxx"):
      pind.append(j)

    if (L[j]=="Lx"):
      aind.append(j)
      print "Lx"
      
  if (len(aind)==0 or len(pind)==0):
    print "Lx or Pxx not found in log file!\n"
    exit(1)

  sf = tstep*float(thermo)
  pind = pind[0]
  aind = aind[0]
  Lxo = dat[0,aind]
 
  strain = nmp.zeros(ntimes,'d')
  # Define the strain(t)
  for i in range(0,ntimes):
    strain[i] = m.log( dat[i,aind] / Lxo )

  otp = open('stress_strain.dat','w')
  for i in range(0,ntimes):
    line = '%f %f \n' % (strain[i], dat[i,pind])
    otp.write(line)
  otp.close()
