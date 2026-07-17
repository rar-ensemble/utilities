import math
import numpy as nmp

def integ(dat):

  N = nmp.shape(dat)[1]
  M = nmp.shape(dat)[0]

  if ( (N==2) and (M>2)):
    dat = dat.T
    N = nmp.shape(dat)[1]
    M = nmp.shape(dat)[0]
  
  de = dat[0][1] - dat[0][0]
  h = (dat[0][N-1] - dat[0][0]) / (N-1)
  
  OM = nmp.zeros(N,'d')

  for i in range(2,N,2):
    OM[i] = OM[i-2] + h * ( dat[1][i-2] + 4 * dat[1][i-1] + dat[1][i]) / 3.0

  for i in range(1,N,2):
    OM[i] = OM[i-1] + h * ( dat[1][i-1] + dat[1][i]) / 2.0

  return OM[N-1]

