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
MAX_STATES = 400
Temperature = 24.054652
HIGHEST_MODE_INDEX = 80



#Initialize all variables
STATES, SADDLE, SADDAT, DATA, Mcon, NE, Nstates, ns, tp, box, ActInd \
    = dimw_initialize.init(MAX_STATES)
nsel = 0


# Check for resume file and continue from previous web
if (os.path.exists('resume.bin')):
  print 'Resuming from file resume.bin!!\n'
  MAX_STATES, STATES, SADDLE, SADDAT, DATA, Mcon, NE, Nstates, ActInd \
      = io.resume()
  print 'Resumed MAX_STATES:', MAX_STATES, 'Ns',Nstates

  nsel, tsel, tsafe, MFPT = ps.pick_state(Mcon, ActInd, NE, Nstates)

  print 'Resuming with nsel:',nsel,'tsel:',tsel,'MFPT:',MFPT,'tsafe:',tsafe
  NE = ps.init_state(nsel, STATES, NE, ActInd, ns, tp, box)


# Main while loop #
while (Nstates < MAX_STATES):

  exmin.explore_min(HIGHEST_MODE_INDEX)

  Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE,  \
      = ana_conn.analyze(nsel,HIGHEST_MODE_INDEX, Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE, Temperature)
  
  io.save_state(MAX_STATES, Nstates, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd)
  
  nsel, tsel, tsafe, MFPT = ps.pick_state(Mcon, ActInd, NE, Nstates)
  
  NE = ps.init_state(nsel, STATES, NE, ActInd, ns, tp, box)
  
  dimw.save_traj(nsel, tsel, tsafe, MFPT, NE, Nstates)
  

