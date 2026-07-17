import math as m
import numpy as nmp
import os
import dimw_utils as dimw
import io_utils as io
import numpy.linalg as lalg
import box0_utils as box0
import solve_master_eqn as meq


T = 24.054

KB = 1.3807E-26 #[kJ / atom / K]
NA = 6.0221415E23 #[ 1 / mol]
RG = NA*KB

KT = KB * T #[Thermal energy in kJ]
RT = RG * T #[Thermal energy in kJ/mol]
Hbar = 1.05457148E-37 #[kJ*s]
H = Hbar * 2 * m.pi #{kJ*s]


MAX, Ns, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd \
    = io.read_state()

x, tp, ns, box = box0.read_all("simul.box0")


DAT2 = nmp.zeros((MAX,4),'d')
SADDAT2 = nmp.zeros((MAX,3),'d')

for i in range(0,MAX):
  DAT2[i,:] = DATA[i,:]
  SADDAT2[i,:] = SADDAT[i,:]


Wsp =nmp.zeros((Ns,3*ns),'d')
Wmn =nmp.zeros((Ns,3*ns),'d')

VSP = nmp.zeros(Ns,'d')
Vmn = nmp.zeros(Ns,'d')

inp = open("Min_eigs.bin","rb")
Wmn = io.read_float_matrix(inp,Ns,3*ns)
inp.close()

inp = open("Sad_eigs.bin","rb")
Wsp = io.read_float_matrix(inp,Ns,3*ns)
inp.close()

for i in range(0,Ns):
  spct = 0
  mnct = 0

  for j in range(3,3*ns):
    if (Wmn[i,j] < 0.0):
      print "Big prob with min", i, Wmn[i,0:4], DATA[i][0]
      name = 'Min_%d_prob.box0' % i
      box0.write_all(ns, STATES[i], tp, box, name)

    #Vmn[i] += m.log( m.sqrt(Wmn[i,j]) )
    Vmn[i] += m.log( m.sqrt(Wmn[i,j]) * 1.0E13 * Hbar / KT )
    mnct += 1
  
  if (i > 0):
    for j in range(4,3*ns):
      if (Wsp[i,j] > 0.001):
        #VSP[i] += m.log( m.sqrt(Wsp[i,j]) )
        VSP[i] += m.log( 1.0E13 * Hbar * m.sqrt(Wsp[i,j]) / KT)
        spct += 1
  if ( (mnct-spct) > 1):
    print "funky min:", mnct, spct, Wsp[i,0:6]

VSP *= RT
Vmn *= RT


Ebar = nmp.zeros((MAX, MAX), 'd')
PEbar = nmp.zeros((MAX, MAX), 'd')
totalct = negct = 0

for i in range(1,Ns):
  totalct += 2
  i1 = int( SADDAT[i,1] )
  i2 = int( SADDAT[i,2] )
  Esp = SADDAT[i,0]
  E1 = DATA[i1,0]
  E2 = DATA[i2,0]

  G1 = E1 + Vmn[i1]
  G2 = E2 + Vmn[i2]
  Gsp = Esp + VSP[i]
  
  print (Gsp-G1),(Gsp-G2),(Esp-E1),(Esp-E2),(VSP[i]-Vmn[i1]), (VSP[i] - Vmn[i2])
  #print Gsp, G1, G2, Esp, E1, E2, VSP[i], Vmn[i1], Vmn[i2]

  Ebar[i1,i2] = Gsp - G1
  Ebar[i2,i1] = Gsp - G2

  DAT2[i2,0] = G2
  DAT2[i1,0] = G1
  SADDAT2[i,0] = Gsp
  
  #print "(",i1, ",",i2,")", Gsp, G1, VSP[i], Vmn[i1], Vmn[i2]

  PEbar[i1,i2] = Esp - E1
  PEbar[i2,i1] = Esp - E2

  if ( (Gsp - G1) < 0):
    negct += 1

  if ( (Gsp - G2) < 0):
    negct += 1

print "Total transitions:", totalct, "Negative", negct

print "Input rates\n", Mcon[0:5, 0:5]

Mcon2 = nmp.zeros((MAX,MAX),'d')
for i in range(0,Ns):
  for j in range(0,Ns):
    if (PEbar[i,j] != 0):
      Mcon2[i,j] = 1.0E-13 * KT * m.exp( -Ebar[i,j] / RT ) / H
    if (PEbar[i,j] != 0):
      Mcon[i,j] =  1.0E-13 * KT * m.exp(-PEbar[i,j] / RT ) / H


print "Calculated rates\n", Mcon[0:5, 0:5]

print "Free energy rates\n", Mcon2[0:5, 0:5]

W2 = meq.form_W(Mcon2,ActInd, NE, Ns)


Q0 = nmp.zeros(NE,'d')
Q0[0] = 1.0
Lam2, S2 = lalg.eig(W2)
Sinv2 = lalg.inv(S2)

W = meq.form_W(Mcon,ActInd, NE, Ns)
Lam, S = lalg.eig(W)
Sinv = lalg.inv(S)

print Lam2

print Lam

MF1 = dimw.calc_MFPT(Q0,S2,Sinv2,Lam2)
MF2 = dimw.calc_MFPT(Q0,S,Sinv,Lam)

print "MFPT:", MF1,MF2

print "safe", dimw.calc_tsafe(Q0,S2,Sinv2,Lam2,MF1), dimw.calc_tsafe(Q0,S,Sinv,Lam,MF2)


io.save_state(MAX,Ns,ns,DAT2,STATES,SADDLE,SADDAT2, Mcon2, NE, ActInd)


