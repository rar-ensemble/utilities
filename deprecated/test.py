import os
import time
import numpy as nmp
import numpy.linalg as lalg
import random as rand
import math
import dimw_utils as dimw
import integrate


W = nmp.zeros((4,4),'d')
Q0 = nmp.zeros((4,1),'d')
Q = nmp.zeros((4,1),'d')
ELam = nmp.zeros((4,4),'d')

#Initial prob. dist
Q0[0] = 0
Q0[1] = 0
Q0[2] = 1
Q0[3] = 0

Q0 = Q0.T

Q0[0][2] = 1.0

# Initialize W
for i in range(0,4):
  W[i,i] = -1.0

W[0,1] = W[0,2] = W[0,3] = 0.25
W[1,0] = W[1,2] = 1.0/3.0
W[2,0] = W[2,1] = W[2,3] = 1.0/3.0
W[3,0] = W[3,2] = 1.0/3.0


# Diagonalize W
Lam, S = lalg.eig(W)

# Invert S
Sinv = lalg.inv(S)

MFPT = dimw.calc_MFPT(Q0,S,Sinv,Lam)

tsafe = dimw.calc_tsafe(Q0,S,Sinv,Lam)

print 'tsafe:',tsafe,'MFPT:', MFPT

