import io_utils as io
import math
import pick_state as ps
import dimw_utils as dimw
import solve_master_eqn as meq
import numpy as nmp
import numpy.linalg as lalg
import random as rand

NP = 5000


MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd \
    = io.read_state()

print 'ns:', ns, 'NE:', NE, 'Nstates:', Ns

W = meq.form_W(Mcon, ActInd, NE, Ns)

Lam, S = lalg.eig(W)
Sinv = lalg.inv(S)

Q0 = nmp.zeros(NE,'d')
Q0[0] = 1.0

dt = 0.025

dat = nmp.zeros((NP,3),'d')
print DATA[0][0]

for i in range(0,NP):
  dat[i][0] = dt*i
  Q = dimw.eval_q(Q0,S,Sinv,Lam,dat[i][0])
  norm = nmp.sum(Q)

  Eavg = 0.0
  for j in range(0,Ns):
    if (ActInd[j] >= 0):
      i1 = ActInd[j]
      Eavg += Q[i1] * DATA[j][0] / norm

  dat[i][1] = Eavg
  dat[i][2] = Q[0] / norm;

print 'Q at tsafe:\n', Q[0]

otp = open('energy_time.dat','w')
for i in range(0,NP):
  line = '%f %f %e\n' % (dat[i][0], dat[i][1], dat[i][2])
  otp.write(line)

otp.close()

