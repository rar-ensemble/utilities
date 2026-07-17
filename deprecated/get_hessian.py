import io_utils as io
import dimw_utils as dimw
import numpy as nmp
import numpy.linalg as lalg
import box0_utils as box0
import os


def calc_hess(ns, x, tp, box):
  os.system("rm hess/*")
  os.system("cp simul.par* ./hess/.")
  os.system("cp simul.input ./hess/.")
  
  box0.write_all(ns,x,tp,box,"./hess/simul.box0")

  os.chdir("./hess")

  os.system("../viklas")
  
  os.chdir("../")



MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd \
    = io.read_state()

print 'ns:', ns, 'NE:', NE, 'Nstates:', Ns

x, tp, ns, box = box0.read_all("simul.box0")

Wsa = nmp.zeros((Ns,3*ns),'d')


for i in range(0,Ns):
  calc_hess(ns,STATES[i],tp,box)

  inp = open("./hess/hessian.output0","rb")
  H = io.read_float_matrix(inp,3*ns,3*ns)
  inp.close()    
  
  W = lalg.eigvalsh(H)
  
  W = nmp.sort(W)
  
  Wsa[i] = W

  if (i%5==0):
    otp = open("Min_eigs.bin","wb")
    io.write_float_matrix(otp,Wsa,Ns,3*ns)
    otp.close()
    print "Checkpoint", i 

  
otp = open("Min_eigs.bin","wb")
io.write_float_matrix(otp,Wsa,Ns,3*ns)
otp.close()

print "Done and done!"



for i in range(1,Ns):
  calc_hess(ns,SADDLE[i],tp,box)

  inp = open("./hess/hessian.output0","rb")
  H = io.read_float_matrix(inp,3*ns,3*ns)
  inp.close()    
  
  W = lalg.eigvalsh(H)
  
  W = nmp.sort(W)
  
  Wsa[i] = W

  if (i%5==0):
    otp = open("Sad_eigs.bin","wb")
    io.write_float_matrix(otp,Wsa,Ns,3*ns)
    otp.close()
    print "Checkpoint", i 

  
otp = open("Sad_eigs.bin","wb")
io.write_float_matrix(otp,Wsa,Ns,3*ns)
otp.close()
print "Done and done!"
