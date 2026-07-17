import numpy as nmp

def remove_jumps(STATES, Mcon, box):
  boxh = box / 2.0

  sh = nmp.shape(STATES)
  Ns = sh[0]
  ns = sh[1]

  for i in range(0,Ns-1):
    for j in range(i+1, Ns):
      
      if (Mcon[i][j] == 0.0):
        continue
  
      dr = STATES[i] - STATES[j]

      for k in range(0,ns):
        if (dr[k][0] > boxh[0]): dr[k][0] -= box[0]
        elif (dr[k][0] < -boxh[0]): dr[k][0] += box[0]

        if (dr[k][1] > boxh[1]): dr[k][1] -= box[1]
        elif (dr[k][1] < -boxh[1]): dr[k][1] += box[1]

        if (dr[k][2] > boxh[2]): dr[k][2] -= box[2]
        elif (dr[k][2] < -boxh[2]): dr[k][2] += box[2]

      STATES[j] = STATES[i] - dr

  return STATES

