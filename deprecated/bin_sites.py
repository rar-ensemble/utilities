import math as m
import numpy as nmp
import pbc_utils as pbc

def bin_sites(x, bx, nbins):
  nx = nbins[0]
  ny = nbins[1]
  nz = nbins[2]

  ns = nmp.shape(x)[0]
  bxh = 0.5*bx
  
  # Bin size in each direction
  drb = bx / nbins
  
  print "Bin sizes:", drb

  binct = nmp.zeros((nx,ny,nz),'d')


  # Put all the particles in the box
  for i in range(0,ns):
    x[i] = pbc.pbc_pos_inbox(x[i], bx)
    
    cbx = int(x[i,0] / drb[0])
    cby = int(x[i,1] / drb[1])
    cbz = int(x[i,2] / drb[2])

    binct[cbx,cby,cbz] += 1.0

  # Define the center of each cube
  binmid = nmp.zeros((nx,ny,nz,3),'d')
  for i in range(0,nx):
    for j in range(0,ny):
      for k in range(0,nz):
        binmid[i,j,k,0] = (i+0.5)*drb[0]
        binmid[i,j,k,1] = (j+0.5)*drb[1]
        binmid[i,j,k,2] = (k+0.5)*drb[2]

  return binct, binmid

    
