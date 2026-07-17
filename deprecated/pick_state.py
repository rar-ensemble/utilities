import os
import box0_utils as box0
import numpy as nmp
import numpy.linalg as lalg
import math
import integrate
import random as rand
import solve_master_eqn as meq
import dimw_utils as dimw

def calc_tselect(Q0, S, Sinv, Lam, Mcon, ActInd, Ns, tsafe):
  dat = nmp.zeros((20000,2),'d')
  dt = tsafe

  maxi = 0


  # This loop calculates f(G,t) using Eq. 4
  for i in range(0,20000):
    t = i*dt
    Q = dimw.eval_q(Q0,S,Sinv,Lam,t)
    dat[i][0] = t
    dat[i][1] = dimw.calc_f(Q,Mcon, ActInd, Ns, t)

    maxi = i

    if (i>10 and dat[i][1] < 1.0e-5):
      break

  dat = dat[0:maxi]
  
  # Normalizes f(G,t) to determine tselect as in step 3 in their paper
  sum = integrate.integ(dat)

  for i in range(0, maxi):
    dat[i][1] *= (1.0 / sum)
 

  # Randomly samples normalized f(G,t) to determine tselect
  tsel = -1
  accept = 0
  delt = dat[1][0] - dat[0][0]

  while (accept==0):
    rx = rand.randrange(0,maxi-1,1)
    ry = rand.random()
    if (ry < dat[rx][1]):
      tsel = dat[rx][0]
      accept = 1
      break


  return tsel
 
def calc_nselect(Q, Mcon, ActInd, Ns, tsel):
  nsel = 0

  f = dimw.calc_f(Q, Mcon, ActInd, Ns, tsel)
  accept = 0

  ry = rand.random()

  P = 0.0

  for i in range(0,Ns):
    if (ActInd[i] >= 0):
      continue

    P += ( dimw.calc_fnG(Q,Mcon,ActInd,Ns,tsel,i) / f)

    if (ry < P):
      nsel = i
      break

#  while (accept==0):
#    # Pick a state
#    rx = rand.randrange(0,Ns,1)
#    
#    # Ensure it's a boundary state
#    if (ActInd[rx] >= 0):
#      continue
#
#    # Accept/reject
#    ry = rand.random()
#    fnG = dimw.calc_fnG(Q, Mcon, ActInd, Ns, tsel, rx)
#    
#    if (ry < fnG / f):
#      nsel = rx
#      accept = 1
#      break

  return nsel

# Calculates tsafe, MFPT, tselect, nselect
def pick_state(Mcon, ActInd, NE, Ns):
  W = meq.form_W(Mcon,ActInd,NE,Ns)
  
  Q0 = nmp.zeros(NE,'d')
  Q0[0] = 1.0
  Lam, S = lalg.eig(W)
  Sinv = lalg.inv(S)
  
  MFPT = dimw.calc_MFPT(Q0,S,Sinv,Lam)
  tsafe = dimw.calc_tsafe(Q0,S,Sinv,Lam,MFPT/10)

  while (tsafe < 0 or tsafe > MFPT):
    tsafe = dimw.calc_tsafe(Q0,S,Sinv,Lam,MFPT*rand.random()/10.0)

  tsel = calc_tselect(Q0, S, Sinv, Lam, Mcon, ActInd, Ns, tsafe)

  Q = dimw.eval_q(Q0,S,Sinv,Lam,tsel)
  
  nsel = calc_nselect(Q, Mcon, ActInd, Ns, tsel)

  return nsel, tsel, tsafe, MFPT


# Initializes the chosen state to be explored
def init_state(nsel, STATES, NE, ActInd, ns, tp, box):
  box0.write_all(ns, STATES[nsel], tp, box, "simul.box0")
  ActInd[nsel] = NE
  NE += 1
  os.system('bash reset.sh')

  return NE
