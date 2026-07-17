import math
import dimw_utils as dimw
import math
import os
import numpy as nmp
import numpy.linalg as lalg
import box0_utils as rb0

# Allocates memory for global variables
# Ensures input files are each present
def init(max_states):
  
  x, tp, ns, box = rb0.read_all('initial.box0')

  SADDLE = nmp.zeros((max_states,ns,3),'d')
  SADDAT = nmp.zeros((max_states,4),'d')
  STATES = nmp.zeros((max_states,ns,3),'d')
  DATA = nmp.zeros((max_states,5),'d')
  Mcon = nmp.zeros((max_states,max_states),'d')
  
  ActInd = nmp.zeros(max_states,'i')
  ActInd[1:max_states] = -1
  ActInd[0] = 0
  
  NE = 1
  Nstates = 1

  os.system('cp initial.box0 simul.box0')

  if (os.path.exists('simul.paramall0')==0 or os.path.exists('simul.paramtype')==0 \
      or os.path.exists('find.input')==0 or os.path.exists('follow.input')==0 \
      or os.path.exists('find-modes')==0 or os.path.exists('follow-modes')==0 ):
    print "Input files not found!!"
    exit(1)
  
  for i in range(0,ns):
    STATES[0][i] = x[i]

  inp = open('initial_e.dat','r')
  ln = inp.readline()
  DATA[0][0] = float(ln)
  inp.close()
  
  return STATES, SADDLE, SADDAT, DATA, Mcon, NE, Nstates, ns, tp, box, ActInd
