import numpy as nmp
import math as m
import pbc_utils as pbc

def g_of_r(x1, x2, box, drbin):
  boxh = 0.5*box
  RMAX = nmp.min(boxh)
  nbins = int(RMAX/drbin) + 2

  gr = nmp.zeros((nbins,2),'d')

  ns1 = nmp.shape(x1)[0]
  ns2 = nmp.shape(x2)[0]

  #  Looping over the same sets of coordinates
  if (ns1==ns2):
    ns1 = nmp.shape(x1)[0]
    ns2 = ns1

    for i in range(0,ns1):
      for j in range(0,ns1):
        if (i==j):
          continue
        mdr = pbc.pbc_mdr(x1[i],x1[j],box,boxh)
        bin = int(mdr / drbin)
        if (bin >= nbins or bin < 0):
          continue

        gr[bin][1] += 1.0

  # Two difference sets of particles
  else:
    ns1 = nmp.shape(x1)[0]
    ns2 = nmp.shape(x2)[0]
    for i in range(0,ns1):
      for j in range(0,ns2):
        mdr = pbc.pbc_mdr(x1[i],x2[j],box,boxh)
        bin = int(mdr / drbin)
        if (bin >= nbins or bin < 0):
          continue

        gr[bin][1] += 1.0


  # Normalize g of r
  tvol = box[0]*box[1]*box[2]
  for i in range(0,nbins):
    dr = (i+0.5)*drbin
    gr[i,0] = dr
    
    dV2 = i*i*i
    dV1 = (i+1)*(i+1)*(i+1)

    norm = 4.0 * drbin**3 * (dV1 - dV2) * m.pi / 3.0
    gr[i,1] *= (tvol / norm / ns1 / ns2)

  return gr


  

