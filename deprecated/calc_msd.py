import numpy as nmp
import dimw_utils as dimw

def msd(POS, box):

  sh = nmp.shape(POS)

  tmax = sh[0]
  ns = sh[1]
  if (sh[2] != 3):
    print 'POS has wrong dimensions in msd!\n'
    exit(1)

  boxh = box / 2.0


  # First remove any jumps across the PBCs
  for t in range(1,tmax):
    dr = POS[t] - POS[t-1]

    for i in range(0,ns):

      if (dr[i][0] > boxh[0]): dr[i][0] -= box[0]
      elif (dr[i][0] < -boxh[0]): dr[i][0] += box[0]

      if (dr[i][1] > boxh[1]): dr[i][1] -= box[1]
      elif (dr[i][1] < -boxh[1]): dr[i][1] += box[1]

      if (dr[i][2] > boxh[2]): dr[i][2] -= box[2]
      elif (dr[i][2] < -boxh[2]): dr[i][2] += box[2]

  msd = nmp.zeros(tmax,'d')

  for delt in range(0, tmax):
    for t in range(0, tmax-delt):
      dr = POS[t+delt] - POS[t]
      dr2 = dr*dr
      msd[delt] = nmp.sum(dr2)



