import numpy as nmp
import os
import io_utils as io
import analyze_connected as ana_conn

MAX_STATES, STATES, SADDLE, SADDAT, DATA, Mcon, NE, Nstates, ActInd \
          = io.resume()




Ns , DATA, STATES, SADDLE, SADDAT, Mcon, NE \
    = ana_conn.analyze(22, 100, 48, DATA, STATES, SADDLE, SADDAT, Mcon, NE, 0.2)
