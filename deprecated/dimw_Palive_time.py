import io_utils as io
import math
import pick_state as ps
import dimw_utils as dimw
import solve_master_eqn as meq
import numpy as nmp
import numpy.linalg as lalg
import random as rand

NP = 200


MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd \
    = io.read_state()

print 'ns:', ns, 'NE:', NE, 'Nstates:', Ns

#Mcon = Mcon / 3.9903

ActInd[200:400] = -1

W = meq.form_W(Mcon, ActInd, NE, Ns)

Lam, S = lalg.eig(W)
Sinv = lalg.inv(S)

Q0 = nmp.zeros(NE,'d')
Q0[0] = 1.0

MFPT = dimw.calc_MFPT(Q0,S,Sinv,Lam)
tsafe = dimw.calc_tsafe(Q0,S,Sinv,Lam,MFPT/10)
while (tsafe < 0 or tsafe > MFPT):
  tsafe = dimw.calc_tsafe(Q0,S,Sinv,Lam,MFPT*rand.random()/10)

print "tsafe:",tsafe, "MFPT:", MFPT

