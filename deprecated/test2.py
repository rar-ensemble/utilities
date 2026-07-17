import os
import time
import numpy as nmp
import numpy.linalg as lalg
import random as rand
import math
import dimw_utils as dimw
import integrate

Nstate = 35
NE = 4
NG = 1

Gstates = nmp.zeros(NG,'i')
Estates = nmp.zeros(NE,'i')

# Active contains flags defining whether states are in E or G #
Active = nmp.zeros(Nstate,'i')
Active[1:5] = 1
ActInd = nmp.zeros(Nstate,'i')
ActInd[0:Nstate] = -1

ind=0
for i in range(0,Nstate):
  if (Active[i]):
    ActInd[i] = ind
    ind += 1

print 'Active:\n',Active

Estates[0] = 1
Estates[1] = 2 
Estates[2] = 3 
Estates[3] = 4
Gstates[0] = 0

# M contains the entire rate constant matrix (sets E and G)
M = nmp.zeros((Nstate,Nstate),'d')

M[1][2] = M[1][3] = M[1][4] = M[1][0] = 0.25
M[2][1] = M[2][3] = M[2][0] =  M[3][1] = M[3][2] = M[3][4] = 1.0/3.0
M[4][1] = M[4][3] = M[4][0] =  M[0][1] = M[0][2] = M[0][4] = 1.0/3.0


W = nmp.zeros((NE,NE),'d')
i1 = i2 = 0

for i in range(0,Nstate):
  
  if (ActInd[i] >= 0):
    i1 = ActInd[i]

  else:
    continue
  
  W[i1][i1] = -nmp.sum(M[i])

  for j in range(0,Nstate):

    if (i==j):
      continue

    elif (ActInd[j] < 0):
      continue
  
    elif (ActInd[j] >= 0):
      i2 = ActInd[j]
      W[i1][i2] = M[i][j]


print W
