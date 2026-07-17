import os
import numpy as nmp
import pbc_utils as pbc
import math as m
import io_utils as io

def legendre2(x):
  return (0.5 * (x*x*3.0 - 1.0))

def make_bond_list(nch, chl, step):
  nbonds = nch*(chl-step)
  bonds = nmp.zeros((nbonds,2),'i')
  ind = 0
  for i in range(0,nch):
    for j in range(0,chl-step):
      i1 = i*chl + j
      i2 = i1 + step
      bonds[ind,0] = i1
      bonds[ind,1] = i2
      ind += 1
  return bonds


def bond_autocorr(x, bx, bonds, NP):
  frs = nmp.shape(x)[0]
  ns = nmp.shape(x)[1]
  nbonds = nmp.shape(bonds)[0]
  bxh = 0.5*bx


  if (nmp.shape(bx)[0] != frs):
    print "ERROR:  Number of frames in box variable does not match position variable\n"
    return 0



  # Pick log-spaced times to aid in fitting later on
  times = nmp.logspace(0,m.log10(frs), NP+3)
  delts = nmp.zeros(NP+3, 'i')

  for i in range(0, NP+3):
    delts[i] = int(times[i])
  delts = nmp.unique(delts)
  NP = nmp.shape(delts)[0]

  # Allocate memory for autocorrelation fn
  bavg = nmp.zeros((NP,2),'d')

  # Main loop
  for t0 in range(0,NP):
    dt = delts[t0]
    bavg[t0,0] = dt
    print t0,"of",NP,"dt =",dt
    
    # Loop over time origins
    for t in range(0, frs-dt):

      #Loop over bonds
      for i in range(0,nbonds):
        i1 = bonds[i,0]
        i2 = bonds[i,1]

        rb1 = pbc.pbc_vdr(x[t,i1], x[t,i2], bx[t], bxh[t])
        rb2 = pbc.pbc_vdr(x[t+dt,i1], x[t+dt,i2], bx[t+dt], bxh[t+dt])

        dot = nmp.dot(rb1, rb2)

        mg1 = m.sqrt(nmp.dot(rb1,rb1))
        mg2 = m.sqrt(nmp.dot(rb2,rb2))

        dot *= (1.0 / mg1 / mg2)

        bavg[t0,1] += legendre2(dot)

  #Normalize the autocorrelator
  for i in range(0,NP):
    dt = delts[i]
    if ( (frs-dt) > 0):
      bavg[i,1] *= (1.0 / nbonds / (frs-dt))

  return bavg
        
def bond_autocorr_c(x, bx, bonds, NP):
  frs = nmp.shape(x)[0]
  ns = nmp.shape(x)[1]
  nbonds = nmp.shape(bonds)[0]
  bxh = 0.5*bx


  if (nmp.shape(bx)[0] != frs):
    print "ERROR:  Number of frames in box variable does not match position variable\n"
    return 0

  # Pick log-spaced times to aid in fitting later on
  times = nmp.logspace(0,m.log10(frs), NP+3)
  delts = nmp.zeros(NP+3, 'i')

  for i in range(0, NP+3):
    delts[i] = int(times[i])
  delts = nmp.unique(delts)
  NP = nmp.shape(delts)[0]

  otp = open('to_c_bond_corr.bin','wb')
  io.write_c_binary(otp,x,bx,delts)
  io.write_int(otp,nbonds)
  io.write_int_matrix(otp,bonds,nbonds,2)
  otp.close()

  os.system("python_bond_autocorr to_c_bond_corr.bin")

  bcorr = io.read_vec_text("bond_corr.out",2,NP)

  return bcorr


