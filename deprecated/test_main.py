import random as rand
import numpy as nmp
import numpy.linalg as lalg
import os
import array
import explore_min as exmin
import analyze_connected as ana_conn
import dimw_utils as dimw
import dimw_initialize
import box0_utils as box0
import solve_master_eqn as meq
import pick_state as ps
import io_utils as io


# Global, predefined variables
MAX_STATES = 100
Temperature = 1.0
HIGHEST_MODE_INDEX = 0


MAX_STATES, STATES, SADDLE, SADDAT, DATA, Mcon, NE, Nstates, ActInd \
          = io.resume()

Ns = Nstates
NE = 2
ActInd[0] = 0
ActInd[1] = 1
ActInd[2:MAX_STATES] = -1
ns = 1
tsafe = 0


W = meq.form_W(Mcon, ActInd, NE, Ns)
Lam, S = lalg.eig(W)
Sinv = lalg.inv(S)
Q0 = nmp.zeros(NE,'d')
Q0[1] = 1.0


while (NE < 98):
  nsel, tsel, tsafe, MFPT = ps.pick_state(Mcon,ActInd, NE, Nstates)

  print tsafe, nsel

  ActInd[nsel] = NE
  dimw.save_traj(nsel, tsel, tsafe, MFPT, NE, Ns)

  NE += 1
  if (NE%10==0):
    io.save_state(MAX_STATES,Nstates,ns,DATA,STATES,SADDLE,SADDAT,Mcon,NE,ActInd)

#  if (NE==4):
#    t = MFPT
#    Ns = Nstates
#
#    Q0 = nmp.zeros(NE,'d')
#    Q0[0] = 1.0
#    
#    W = meq.form_W(Mcon,ActInd,NE,Ns)
#    Lam, S = lalg.eig(W)
#    Sinv = lalg.inv(S)
#
#    Q = dimw.eval_q(Q0,S,Sinv,Lam,t)
#    print t
#    for i in range(0,NE):
#      if (ActInd[i]<0):
#        print "fnG", i, ":", dimw.calc_fnG(Q, Mcon, ActInd, Ns, t, i)
#        break
#    exit(0)


io.save_state(MAX_STATES,Nstates,ns,DATA,STATES,SADDLE,SADDAT,Mcon,NE,ActInd)

#ns = 1
#
#SADDLE = nmp.zeros((MAX_STATES,ns,3),'d')
#SADDAT = nmp.zeros((MAX_STATES,3),'d')
#STATES = nmp.zeros((MAX_STATES,ns,3),'d')
#DATA = nmp.zeros((MAX_STATES,4),'d')
#
#
#inp = open("test_system.bin","rb")
#Mcon = io.read_float_matrix(inp,100,100)
#
#Ns = 100
#NE = 100
#ActInd = nmp.zeros(100,'d')
#for i in range(0,NE):
#  ActInd[i] = i
#  i1, i2 = dimw.unwrap_inds(i,10,10)
#  
#  E = 0.0001*i1 + 0.001*i2
#
#  DATA[i,0] = E
#  
#io.save_state(MAX_STATES, Ns, 1, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd)




