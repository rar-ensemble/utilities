import math as m
import pbc_utils as pbc
import numpy as nmp

def sort_cluster(W, mx):
  for i in range(0, mx):

    j = 0
    while (j < mx):
      if (i==j):
        j += 1
        continue
      
      flag = 1
      if (W[i,j] == 1):

        sum = 0
        for k in range(0,mx):
          sum += W[j,k]

        if (sum > 0):
          for k in range(0,mx):
            W[i,k] += W[j,k]

            if (W[i,k] > 1):
              W[i,k] = 1

            W[j,k] = 0

          j=0
          flag = 0

      if (flag==1):
        j += 1
  return W


def empty_bins(binct, binmid, p1, bx, rcut):
  nx = nmp.shape(binct)
  
  bh = 0.5*bx

  mtbins = nmp.zeros(nx,'i')

  unw = nmp.reshape(binct,(-1,1))

  N = nmp.shape(unw)[0]

  ns = nmp.shape(p1)[0]

  # First, loop over all of the bins and denote those
  # that are not near other particles
  num_MT = 0
  for i in range(0,nx[0]):
    for j in range(0,nx[1]):
      for k in range(0,nx[2]):
        
        flag = 1
        for p in range(0,ns):
          if (pbc.pbc_mdr(binmid[i,j,k], p1[p], bx, bh) < rcut):
            flag = 0
            break

        if (flag==1):
          mtbins[i,j,k] = 1
          num_MT += 1

  print num_MT, "of", N, "qualifying bins!"
  return mtbins

def volume_clusters(binct, binmid, p1, bx):
  nx = nmp.shape(binct)
  
  bh = 0.5*bx

  mtbins = nmp.zeros(nx,'i')

  unw = nmp.reshape(binct,(-1,1))

  N = nmp.shape(unw)[0]

  ns = nmp.shape(p1)[0]

  # First, loop over all of the bins and denote those
  # that are not near other particles
  num_MT = 0
  for i in range(0,nx[0]):
    for j in range(0,nx[1]):
      for k in range(0,nx[2]):
        
        flag = 1
        for p in range(0,ns):
          if (pbc.pbc_mdr(binmid[i,j,k], p1[p], bx, bh) < 4.0):
            flag = 0
            break

        if (flag==1):
          mtbins[i,j,k] = 1
          num_MT += 1
  
  print num_MT, "of", N, "qualifying bins!"

  inds = nmp.zeros(nx,'i')
  un_inds = nmp.zeros((N,3),'i')
  unMT = nmp.zeros(N,'i')
  ind = 0
  for i in range(0,nx[0]):
    for j in range(0,nx[1]):
      for k in range(0,nx[2]):
        inds[i,j,k] = ind
        un_inds[ind,0] = i
        un_inds[ind,1] = j
        un_inds[ind,2] = k

        unMT[ind] = mtbins[i,j,k]
        ind += 1

  
  clust = nmp.zeros((N,N),'i')

  for i in range(0,nx[0]):
    for j in range(0,nx[1]):
      for k in range(0,nx[2]):
        
        # If this bin is empty, keep going
        if (mtbins[i,j,k] == 0):
          continue


        i1 = inds[i,j,k]

        # +X direction
        ii = i + 1
        if (ii >= nx[0]):
          ii -= nx[0]
        if (mtbins[ii,j,k] == 1):
          i2 = inds[ii,j,k]
          clust[i1,i2] = clust[i2,i1] = 1

        # -X direction
        ii = i - 1
        if (ii < 0):
          ii += nx[0]
        if (mtbins[ii,j,k] == 1):
          i2 = inds[ii,j,k]
          clust[i1,i2] = clust[i2,i1] = 1

        # +Y direction
        ii = j + 1
        if (ii >= nx[1]):
          ii -= nx[1]
        if (mtbins[i,ii,k] == 1):
          i2 = inds[i,ii,k]
          clust[i1,i2] = clust[i2,i1] = 1

        # -Y direction
        ii = j - 1
        if (ii < 0):
          ii += nx[1]
        if (mtbins[i,ii,k] == 1):
          i2 = inds[i,ii,k]
          clust[i1,i2] = clust[i2,i1] = 1

        # +Z direction
        ii = k + 1
        if (ii >= nx[2]):
          ii -= nx[2]
        if (mtbins[i,j,ii] == 1):
          i2 = inds[i,j,ii]
          clust[i1,i2] = clust[i2,i1] = 1

        # -Y direction
        ii = j - 1
        if (ii < 0):
          ii += nx[2]
        if (mtbins[i,j,ii] == 1):
          i2 = inds[i,j,ii]
          clust[i1,i2] = clust[i2,i1] = 1

        
  return clust, mtbins















