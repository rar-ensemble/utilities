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

  Lind = []
  pind = []
  for j in range(0,len(L)):
    
    if (L[j]=="Pxx"):
      pind.append(j)
    
    if (L[j]=="Pyy"):
      pind.append(j)
    
    if (L[j]=="Pzz"):
      pind.append(j)

    if (L[j]=="Lx"):
      Lind.append(j)

    if (L[j]=="Ly"):
      Lind.append(j)

    if (L[j]=="Lz"):
      Lind.append(j)
      
  if (len(Lind)<3 or len(pind)<3):
    print "Lx or Pxx not found in log file!\n"
    exit(1)

  sf = tstep*float(thermo)

  Px = pind[0]
  Py = pind[1]
  Pz = pind[2]

  Lx = Lind[0]
  Ly = Lind[1] 
  Lz = Lind[2]

  Press = (dat[:,Px] + dat[:,Py] + dat[:,Pz]) / 3.0
  Vol = (dat[:,Lx] * dat[:,Ly] * dat[:,Lz])
  
  Vo = Vol[0]

  otp = open('press_del-vol.dat','w')
  for i in range(0,ntimes):
    line = '%f %f \n' % ((Vol[i] / Vo - 1.0), Press[i])
    otp.write(line)
  otp.close()
