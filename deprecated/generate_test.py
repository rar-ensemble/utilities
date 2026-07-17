import random as rand
import numpy as nmp
import math as m
import io_utils as io

Umin = nmp.zeros((10,10),'d')
Uvec = nmp.zeros((10*10,3),'d')

for i in range(0,10):
  for j in range(0,10):
    Umin[i,j] = 0.0001*i + 0.001*j
    Uvec[i*10+j,0] = Umin[i,j]
    Uvec[i*10+j,1] = i
    Uvec[i*10+j,2] = j


SPE = nmp.zeros((500,3),'d')
nsad = 0.0

for i in range(0,10):
  for j in range(0,10):
    i1 = i*10+j
    
    if ( (i+1)<10 ):
      d = 1.0
      
      i2 = (i+1)*10 + j
      
      SPE[nsad,0] = max(Uvec[i1,0], Uvec[i2,0]) + \
          abs( (Uvec[i1,0] - Uvec[i2,0]) ) * (0.001 + 10*rand.random()*d)
      SPE[nsad,1] = i1
      SPE[nsad,2] = i2
      
      nsad += 1.0


    if ( (j+1) < 10):
      d = 2.0
      
      i2 = (i)*10 + j+1
      
      SPE[nsad,0] = max(Uvec[i1,0], Uvec[i2,0]) + \
          abs( (Uvec[i1,0] - Uvec[i2,0]) ) * (0.001 + 10*rand.random()*d)
      SPE[nsad,1] = i1
      SPE[nsad,2] = i2
      
      nsad += 1.0
print "determined", nsad, "saddles"

Kmat = nmp.zeros((100,100),'d')
for i in range(0,nsad):
  i1 = SPE[i,1]
  i2 = SPE[i,2]
  
  print (SPE[i,0] - Uvec[i1,0]), (SPE[i,0] - Uvec[i2,0])
  Kmat[i1,i2] = m.exp(- (SPE[i,0] - Uvec[i1,0]))
  Kmat[i2,i1] = m.exp(- (SPE[i,0] - Uvec[i2,0]))

otp = open("test_system.bin","wb")
io.write_float_matrix(otp,Kmat,100,100)
otp.close()
