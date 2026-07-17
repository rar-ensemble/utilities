import io_utils as io
import numpy as nmp




MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd \
    = io.read_state()

print 'ns:', ns, 'NE:', NE, 'Nstates:', Ns, 'Old Max:', MAX

MAX2 = MAX*2

DAT2 = nmp.zeros((MAX2,5),'d')
DAT2[0:MAX,0:4] = DATA

STATES2 = nmp.zeros((MAX2,ns,3),'d')
STATES2[0:MAX,0:ns,0:3] = STATES

SADDLE2 = nmp.zeros((MAX2,ns,3),'d')
SADDLE2[0:MAX,0:ns,0:3] = SADDLE

SADDAT2 = nmp.zeros((MAX2,4),'d')
SADDAT2[0:MAX,0:3] = SADDAT

M2 = nmp.zeros((MAX2,MAX2),'d')
M2[0:MAX,0:MAX] = Mcon

Act2 = nmp.zeros(MAX2,'d')
Act2[0:MAX] = ActInd
Act2[MAX:MAX2] = -1



io.save_state(MAX2, Ns, ns, DAT2, STATES2, SADDLE2, SADDAT2, M2, NE, Act2)
