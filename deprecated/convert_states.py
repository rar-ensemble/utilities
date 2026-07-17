import io_utils as io


MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT2, Mcon, NE, ActInd \
    = io.read_state()

print 'ns:', ns, 'NE:', NE, 'Nstates:', Ns

Mcon = Mcon / 3.9903

io.save_state(MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT2, Mcon, NE, ActInd)
