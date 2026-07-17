import numpy as nmp
import math as m

def shear_nonaffine( x1 , x2 , bx , sh1 , sh2 ):
  ns = nmp.shape( x1 ) [0]
  dshear = sh2 - sh1 
  
  disp = nmp.zeros( (ns,4) , 'd' )

  for i in range( 0 , ns ):
    for j in range( 0 , 3 ) :
      while ( x1[i,j] > 0.5 * bx[j] ):
        x1[i,j] -= bx[j]

      while ( x1[i,j] < -0.5* bx[j] ):
        x1[i,j] += bx[j]

    x1[i,0] += x1[i,1] * dshear

    for j in range( 0 , 3 ):
      disp[i,j] = x2[i,j] - x1[i,j]

      while ( disp[i,j] > 0.5 * bx[j] ): 
        disp[i,j] -= bx[j]
      
      while ( disp[i,j] < -0.5* bx[j] ):
        disp[i,j] += bx[j]

      disp[i,3] += disp[i,j]*disp[i,j]

  return disp

def pratio(x1, x2, bx1, bx2):
  ns = nmp.shape(x1)[0]
  if (nmp.shape(x2)[0] != ns):
    print "shape of x1 and x2 do not agree!"
    return -1.0

  s1 = nmp.zeros((ns,3),'d')
  s2 = nmp.zeros((ns,3),'d')
  af = nmp.zeros((ns,3),'d')

  ds = nmp.zeros(3,'d')
  cmds = nmp.zeros(3,'d')

  # Calculate scaled coords, shifts in center of mass
  for i in range(0,ns):
    for j in range(0,3):
      s1[i,j] = x1[i,j] / bx1[j]
      s2[i,j] = x2[i,j] / bx2[j]

      ds[j] = s2[i,j] - s1[i,j]

      if (ds[j] > 0.5):
        ds[j] -= 1.0
      elif (ds[j] <= -0.5):
        ds[j] += 1.0

      s2[i,j] = s1[i,j] + ds[j]

      cmds[j] += ds[j]


  cmds *= (1.0 / ns)

  
  # Calculate nonaffine displacements
  pratio = 0.0
  av_dr2 = 0.0
  av_dr4 = 0.0
  for i in range(0,ns):
    ds = s2[i,:] - s1[i,:] - cmds

    dr = ds * bx2

    mdr2 = nmp.dot(dr,dr)
    
    av_dr2 += mdr2
    av_dr4 += (mdr2*mdr2)

  pratio = av_dr2*av_dr2 / (ns * av_dr4)

  return pratio


def nonaff(x1, x2, bx1, bx2):
  ns = nmp.shape(x1)[0]
  if (nmp.shape(x2)[0] != ns):
    print "shape of x1 and x2 do not agree!"
    return -1.0

  s1 = nmp.zeros((ns,3),'d')
  s2 = nmp.zeros((ns,3),'d')
  af = nmp.zeros((ns,3),'d')

  ds = nmp.zeros(3,'d')
  cmds = nmp.zeros(3,'d')

  # Calculate scaled coords, shifts in center of mass
  for i in range(0,ns):
    for j in range(0,3):
      s1[i,j] = x1[i,j] / bx1[j]
      s2[i,j] = x2[i,j] / bx2[j]

      ds[j] = s2[i,j] - s1[i,j]

      if (ds[j] > 0.5):
        ds[j] -= 1.0
      elif (ds[j] <= -0.5):
        ds[j] += 1.0

      s2[i,j] = s1[i,j] + ds[j]

      cmds[j] += ds[j]


    cmds *= (1.0 / ns)

  
  # Calculate nonaffine displacements
  pratio = 0.0
  av_dr2 = 0.0
  av_dr4 = 0.0
  for i in range(0,ns):
    ds = s2[i,:] - s1[i,:] 
    if (ns > 1):
      ds -= cmds

    dr = ds * bx2

    mdr2 = nmp.dot(dr,dr)
    
    av_dr2 += mdr2
    av_dr4 += (mdr2*mdr2)

  return av_dr2/ns



    

    


    

def nonaff_noisy(x1, x2, bx1, bx2):
  ns = nmp.shape(x1)[0]
  if (nmp.shape(x2)[0] != ns):
    print "shape of x1 and x2 do not agree!"
    return -1.0

  s1 = nmp.zeros((ns,3),'d')
  s2 = nmp.zeros((ns,3),'d')
  af = nmp.zeros((ns,3),'d')

  ds = nmp.zeros(3,'d')
  cmds = nmp.zeros(3,'d')

  # Calculate scaled coords, shifts in center of mass
  for i in range(0,ns):
    for j in range(0,3):
      s1[i,j] = x1[i,j] / bx1[j]
      s2[i,j] = x2[i,j] / bx2[j]

      ds[j] = s2[i,j] - s1[i,j]

      if (ds[j] > 0.5):
        ds[j] -= 1.0
      elif (ds[j] <= -0.5):
        ds[j] += 1.0

      s2[i,j] = s1[i,j] + ds[j]

      cmds[j] += ds[j]


    cmds *= (1.0 / ns)

  
  # Calculate nonaffine displacements
  pratio = 0.0
  av_dr2 = 0.0
  av_dr4 = 0.0
  for i in range(0,ns):
    ds = s2[i,:] - s1[i,:] 
    if (ns > 1):
      ds -= cmds

    dr = ds * bx2

    mdr2 = nmp.dot(dr,dr)
    
    print i, mdr2

    av_dr2 += mdr2
    av_dr4 += (mdr2*mdr2)

  return av_dr2/ns



    

    
