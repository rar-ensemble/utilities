import os
import time
import integrate
import numpy as nmp
import numpy.linalg as lalg
import random as rand
import math
import dimw_utils as dimw
import box0_utils as box0
import io_utils as io

dt = 0.05
Etol = 0.001
R2tol = 0.01

Planck = 6.626068E-24 #Includes factor of 10^13 for LJ units, 10^-3 for kJ
Boltz = 1.3806503E-26 #kJ/K
NA = 6.0221415E23

def eval_q(Q0,S,Sinv,Lam,t):

  N = nmp.shape(S)[0]
  ELam = nmp.zeros((N,N),'d')

  for i in range(0,N):
    ELam[i,i] = math.exp( (Lam[i]*t) )

  Q = nmp.dot(Q0,nmp.dot(S,nmp.dot(ELam,Sinv)))

  return Q

def eval_q2(Q0,R,L,Lam,t):

  N = nmp.shape(S)[0]
  ELam = nmp.zeros((N,N),'d')

  for i in range(0,N):
    ELam[i,i] = math.exp( (Lam[i]*t) )

  #Q = nmp.dot(Q0,nmp.dot(S,nmp.dot(ELam,Sinv)))
  Q = nmp.dot(Q0, nmp.dot(R, nmp.dot(ELam, L)))

  return Q


def Palive(Q):

  return nmp.sum(Q)

def calc_tsafe(Q0, S, Sinv, Lam, MFPT):
  t = MFPT
  Palive = 0.0
  while (Palive < 0.999):
    t *= 0.98
    Palive = dimw.Palive(dimw.eval_q(Q0,S,Sinv,Lam,t))

  return t


#def calc_tsafe(Q0,S,Sinv,Lam,guess):
#  inft = dt / 100.0
#    
#  
#  if (guess < 0):
#    print 'Your guess for tsafe sucks and I quit!\n'
#    exit(1)
#
#
#  tnew = guess
#  
#  err = 500
#
#  ct = 0
#
#  while (err > 0.00005) and ct < 20:
#    
#    if (tnew < 0.0):
#      tnew = guess*rand.random()
#
#    Q = eval_q(Q0,S,Sinv,Lam,tnew)
#    Pnew = Palive(Q)
#    Errnew = Pnew - 0.999    
#    
#    delt = tnew + inft
#    Q = eval_q(Q0,S,Sinv,Lam,delt)
#    #Pf = (Palive(Q) - 0.999)
#    Pf = math.fabs(Palive(Q) - 0.999)
#
#    delt = tnew - inft
#    Q = eval_q(Q0,S,Sinv,Lam,delt)
#    #Pb = (Palive(Q) - 0.999)
#    Pb = math.fabs(Palive(Q) - 0.999)
#
#    dP = (Pf - Pb) / (2.0*inft)
#    
#    told = tnew 
#
#    # tchange checks to see how big of a change in t is occurring
#    # if too large, scales tchange to smaller value
#    tchange = Errnew / dP
#
#    if ( math.fabs(tchange / told) > 3.0):
#      if (tchange < 0):
#        tchange = -3.0*told
#      else:
#        tchange = 3.0*told
#
#    tnew = told + tchange
#    
#    err = math.fabs(Errnew)
#
#    ct += 1
#
#    if (ct==2500):
#      tnew = -1.0
#      print "I suck and I'm starting over!\n"
#      break
#
#  return tnew


def calc_MFPT2(Q0,R,L,Lam):
  
  N = 1000
  dat = nmp.zeros((2,N),'d')
  dat[1,N-1] = 1.0

  delt = dt

  while (dat[1,N-1] > 1.0E-4):
    for i in range(0,N):
      t = delt*float(i)
      
      dat[0][i] = t
  
      Q = eval_q2(Q0,R,L,Lam,t)
  
      dat[1][i] = Palive(Q)
    delt *= 10

  return integrate.integ(dat)


def calc_MFPT(Q0,S,Sinv,Lam):
  
  N = 1000
  dat = nmp.zeros((2,N),'d')
  dat[1,N-1] = 1.0

  delt = dt

  while (dat[1,N-1] > 1.0E-4):
    for i in range(0,N):
      t = delt*float(i)
      
      dat[0][i] = t
  
      Q = eval_q(Q0,S,Sinv,Lam,t)
  
      dat[1][i] = Palive(Q)
    delt *= 10

  return integrate.integ(dat)

# Calculates f_nG(G,t) for a given distribution
# This quantity is the rate into boundary state 
# nG from the active states
def calc_fnG(Q, Mcon, ActInd, Ns, t, Gind):
  if (ActInd[Gind] >= 0):
    print 'ERROR!!  cannot calculate fnGt for an active state!\n'
    exit(1)


  fnG = 0.0
  i1 = -1
  for i in range(0, Ns):
    if (ActInd[i] < 0):
      continue
    
    i1 = ActInd[i]  # Need the 'active' index so we know which Q to take

    fnG += Q[i1] * Mcon[i][Gind]

    
  return fnG


# Sums over each boundary state to find total
# Flux into the boundary
def calc_f(Q, Mcon, ActInd, Ns, t):
  f = 0.0
  for i in range(0,Ns):
    if (ActInd[i] < 0):
      f += calc_fnG(Q, Mcon, ActInd, Ns, t, i)
  
  return f



# box0_rmsd - returns the mean squared between two given
# configurations
def box0_rmsd(name1, name2):
  # Run rmsd code on first min
  cmd = "box0_rmsd %s %s >py-stdout" % (name1, name2)
  os.system(cmd)
  
  inp = open('py-stdout','r')
  L = inp.readline()
  L = inp.readline()
  L = inp.readline()
  L = L.split()
  r1m1 = float(L[1])

  inp.close()

  return r1m1



# valid_min checks energies in mode_0_trans_data.dat
# to see if two distinct minima were found, and one of 
# them is the initial minimum.  
# returns flag, E1, E2, Es
# flag = flag for validity of min.
# E1 = initial minimum energy
# E2 = second minimum energy
# Es = saddle energy
def valid_min(i):
  
  name = '%d/mode_0_trans_data.dat' % i
  inp = open(name,'r')

  # Read initial energy
  L = inp.readline()
  L = L.split()
  Ei = float(L[1])

  # Read min. 1 energy
  L = inp.readline()
  L = L.split()
  E1 = float(L[1])

  # Read saddle energy
  L = inp.readline()
  L = L.split()
  Es = float(L[1])

  # Read min. 2 energy
  L = inp.readline()
  L = L.split()
  E2 = float(L[1])

  inp.close()


  # These two flags detect whether this saddle
  # a) Found two distinct minima
  # b) Returned to the initial minimum
  # c) If swap==1, the file index 2 corresponds to the starting minimum
  two_min = 0
  old_min = 0
  flag = 0
  swap = 0

  if ( math.fabs( (E1-E2) ) > Etol ):
    two_min = 1

  if ( math.fabs( (E1-Ei) ) < Etol ):
    old_min = 1
  else:
    E1, E2 = E2, E1
    swap = 1
    if ( math.fabs( (E1-Ei) ) < Etol ):
      old_min = 1
   
  # If the configurations have passed the energy tests,
  # Check the r^2 to ensure one config. went back to starting min
  if (old_min):
    
    n1 = '%d/simul.box0' % i
    if (swap==1):
      n2 = '%d/mode_0_min_%d.box0.sav' % (i,2)
    else:
      n2 = '%d/mode_0_min_%d.box0.sav' % (i,1)
    r1 = box0_rmsd(n1,n2)

    if (r1 > R2tol):
      old_min = 0
   
  if (two_min and old_min):
    flag = 1
    
  return flag, swap, E1, E2, Es

 
# min_rmsd - returns the first, second, and third
# moments of the displacement of a new minimum from
# the initial configuration
def min_rmsd(i, swap):

  # Relevant file names
  m1name = '%d/mode_0_min_1.box0.sav' % i
  m2name = '%d/mode_0_min_2.box0.sav' % i

  if (swap==1):
    m2name = '%d/mode_0_min_1.box0.sav' % i
    m1name = '%d/mode_0_min_2.box0.sav' % i

  iname = 'initial.box0'
  
#  # Run rmsd code on first min
#  cmd = "box0_rmsd %s %s >py-stdout" % (iname,m1name)
#  os.system(cmd)
#  
#  # Read in data from first min
#  inp = open('py-stdout','r')
#  L = inp.readline()
#  L = inp.readline()
#  L = inp.readline()
#  L = L.split()
#  r1m1 = float(L[1])
#  L = inp.readline()
#  L = L.split()
#  r1m2 = float(L[1])
#  L = inp.readline()
#  L = L.split()
#  r1m3 = float(L[1])
#  inp.close()


  # Run rmsd code on second min
  cmd = "box0_rmsd %s %s >py-stdout" % (iname,m2name)
  os.system(cmd)
  
  # Read in data from second min
  inp = open('py-stdout','r')
  L = inp.readline()
  L = inp.readline()
  L = inp.readline()
  L = L.split()
  r2m1 = float(L[1])
  L = inp.readline()
  L = L.split()
  r2m2 = float(L[1])
  L = inp.readline()
  L = L.split()
  r2m3 = float(L[1])
  inp.close()


  return r2m1,r2m2,r2m3

def new_min(E2, Es, r1, r2, r3, Nstates, DATA, SADDAT):
  
  flag = -1
  for i in range(0,Nstates):
    dE1 = math.fabs(DATA[i][0] - E2)
    if (i > 0):
      dE2 = math.fabs(SADDAT[i][0] - Es)
    else:
      dE2 = 2342.0
    dr1 = math.fabs(DATA[i][1] - r1)
    dr2 = math.fabs(DATA[i][2] - r2)
    dr3 = math.fabs(DATA[i][3] - r3)
    
    ct = 0

    if (dE1 < Etol):
      ct += 1
    if (dE2 < Etol):
      ct += 1
    if (dr1 < R2tol):
      ct += 1
    if (dr2 < R2tol):
      ct += 1
    if (dr3 < R2tol):
      ct += 1
    
    if (ct > 2):
      flag = i
      break

    #if (dE1 < Etol and dr1 < R2tol and dr2 < R2tol and dr3 < R2tol):
    #  flag = i
    #  break
      
    
  return flag

def add_connection(curi, conn, swap, Mcon, Es, E2, E1, Temp):
  kT = Boltz*Temp
  RT = kT*NA
  rate = kT * math.exp( -(Es-E1) / RT) / Planck
  Mcon[curi][conn] = rate

  rate = kT * math.exp( -(Es-E2) / RT) / Planck
  Mcon[conn][curi] = rate
 
  return Mcon

def get_vibrations(cur_state, ns, kT):
  Vm1 = 0.0
  Vm2 = 0.0
  Vsp = 0.0

  name = '%d/min_1_hessian.bin' % cur_state
  inp = open(name,'rb')
  H = io.read_float_matrix(inp,3*ns,3*ns)
  Wm1 = lalg.eigvalsh(H)
  inp.close()

  name = '%d/saddle_hessian.bin' % cur_state
  inp = open(name,'rb')
  H = io.read_float_matrix(inp,3*ns,3*ns)
  Wsp = lalg.eigvalsh(H)
  inp.close()

  name = '%d/min_2_hessian.bin' % cur_state
  inp = open(name,'rb')
  H = io.read_float_matrix(inp,3*ns,3*ns)
  Wm2 = lalg.eigvalsh(H)
  inp.close()

  for i in range(3,3*ns):
    Vm1 += math.log( Planck * math.sqrt(Wm1[i]) / kT)
    Vm2 += math.log( Planck * math.sqrt(Wm2[i]) / kT)

    if ( (i >= 4) and (Wsp[i] > 0.001) ):
      Vsp += math.log( Planck * math.sqrt(Wsp[i]) / kT)

  RT = kT*NA
  Vsp *= RT
  Vm1 *= RT
  Vm2 *= RT

  return [Vm1, Vm2, Vsp]

# Called from analyze, adds a newly discovered state to the boundary states
def add_state(cur_state,dir, swap, Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE, dat):
  E1 = dat[0]
  E2 = dat[1]
  Es = dat[2]
  r1 = dat[3]
  r2 = dat[4]
  r3 = dat[5]
  Temp = dat[6]
  
  kT = Boltz*Temp
  RT = kT*NA
  
  Ns = Nstates
  
  DATA[Ns][0] = E2
  DATA[Ns][1] = r1
  DATA[Ns][2] = r2
  DATA[Ns][3] = r3

  name = '%d/mode_0_saddle.box0.sav' % dir
  SADDLE[Ns], tp, nsites, box = box0.read_all(name)
  SADDAT[Ns][0] = Es
  SADDAT[Ns][1] = cur_state
  SADDAT[Ns][2] = Ns

  name = '%d/mode_0_min_2.box0.sav' % dir
  if (swap):
    name = '%d/mode_0_min_1.box0.sav' % dir

  STATES[Ns], tp, nsites, box = box0.read_all(name)


  # Check for Hessian files - if they exist,
  # read them in and calculate vibrational entropy
  nm1 = '%d/min_1_hessian.bin' % dir
  nm2 = '%d/min_2_hessian.bin' % dir
  nsp = '%d/saddle_hessian.bin' % dir
  
  FE = [0.0, 0.0, 0.0]
  if (os.path.exists(nm1) and os.path.exists(nm2) and os.path.exists(nsp)):
    FE = get_vibrations(dir, nsites, kT)

  DATA[Ns,4] = FE[1]
  if (swap):
    DATA[Ns,4] = FE[0]
  SADDAT[Ns,3] = FE[2]

  Es += SADDAT[Ns,3]
  E1 += DATA[cur_state,4]
  E2 += DATA[Ns, 4]
 
  # Rate to leave current state
  rate = kT * math.exp( -(Es-E1) / RT) / Planck
  Mcon[cur_state][Ns] = rate

  # Rate into current state
  rate = kT * math.exp( -(Es-E2) / RT) / Planck
  Mcon[Ns][cur_state] = rate

  Nstates += 1

  return Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE


def save_traj(nsel, tsel, tsafe, MFPT, NE, Nstates):
  name = 'dimw_traj.output'
  
  if (NE==1):
    otp = open(name,'w')
    line = '%d %d %d %f %f %f\n' % (NE,Nstates,nsel,tsel,tsafe,MFPT)
    otp.write(line)
    otp.close()
  else:
    otp = open(name,'a')
    line = '%d %d %d %f %f %f\n' % (NE,Nstates,nsel,tsel,tsafe,MFPT)
    otp.write(line)
    otp.close()







  
