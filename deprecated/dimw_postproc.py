import io_utils as io
import math
import pick_state as ps
import dimw_utils as dimw
import solve_master_eqn as meq
import numpy as nmp
import numpy.linalg as lalg


MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd \
    = io.read_state()

print 'ns:', ns, 'NE:', NE, 'Nstates:', Ns

W = meq.form_W(Mcon, ActInd, NE, Ns)

Lam, S = lalg.eig(W)
Sinv = lalg.inv(S)

Q0 = nmp.zeros(NE,'d')
Q0[0] = 1.0

print 'Lam\n', Lam

print 'safe time:', dimw.calc_tsafe(Q0,S,Sinv,Lam,MFPT/10)

nsel, tsel, tsafe, MFPT = ps.pick_state(Mcon, ActInd, NE, Ns)
print ActInd[0:Ns]
print nsel, tsel, tsafe, MFPT
