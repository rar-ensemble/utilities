import numpy as nmp
import pbc_utils as pbc

def min_dist_hist(x1, x2, bx):
  bh = 0.5*bx

  ns1 = nmp.shape(x1)[0]
  ns2 = nmp.shape(x2)[0]

  min = nmp.zeros(ns1,'d')

  for i in range(0,ns1):
    min[i] = 555.2
    for j in range(0,ns2):
      if (ns1==ns2 and i==j):
        continue
      dr = pbc.pbc_mdr(x1[i],x2[j],bx,bh)
      if (dr < min[i]):
        min[i] = dr

  return min
