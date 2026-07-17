import math
import os
import numpy as nmp
import numpy.linalg as lalg
import random as rand
import dimw_utils as dimw
import integrate

def form_W(Mcon, ActInd, NE, Ns):
  W = nmp.zeros((NE,NE),'d')
 
  i1 = i2 = 0
  for i in range(0,Ns):
    if (ActInd[i] >= 0):
      i1 = ActInd[i]
    else:
      continue
    
    W[i1][i1] = -nmp.sum(Mcon[i])

    for j in range(0,Ns):
      if (i==j):
        continue

      elif (ActInd[j] < 0):
        continue

      else:
        i2 = ActInd[j]
        W[i1][i2] = Mcon[i][j]

  return W
