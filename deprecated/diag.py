import numpy as nmp
import numpy.linalg as lalg
import random as rand

A = nmp.zeros((2112,2112),'d')
for i in range(0,2112):
  for j in range(0,2112):
    A[i,j] = rand.random()

print nmp.shape(A)

Lam, S = lalg.eig(A)
