import os
import math as m
import io_utils as io
import numpy as nmp
import psyco

psyco.full()
def all(x, bx, q, NP):

  frs = nmp.shape(x)[0]
  ns = nmp.shape(x)[1]
  bxh = 0.5*bx


  if (nmp.shape(bx)[0] != frs):
    print "ERROR:  Number of frames in box variable does not match position variable\n"
    return 0



  # Pick log-spaced times to aid in fitting later on
  times = nmp.logspace(0,m.log10(frs), NP+3)
  delts = nmp.zeros(NP+3, 'i')

  for i in range(0, NP+3):
    delts[i] = int(times[i])
  delts = nmp.unique(delts)
  NP = nmp.shape(delts)[0]

  otp = open('to_c_dynamics.bin','wb')
  io.write_c_binary(otp,x,bx,delts)
  io.write_int(otp,q[0])
  io.write_int(otp,q[1])
  io.write_int(otp,q[2])
  otp.close()

  os.system('python_dynamics to_c_dynamics.bin')

  # Read in output
  dat = nmp.zeros((NP,4),'d')
  inp = open('dynamics.out','r')
  for i in range(0,NP):
    line = inp.readline().split()
    dat[i,0] = float(line[0])
    dat[i,1] = float(line[1])
    dat[i,2] = float(line[2])
    dat[i,3] = float(line[3])
  inp.close()

  return dat


def nonaff_all(x, bx, q, NP):

  frs = nmp.shape(x)[0]
  ns = nmp.shape(x)[1]
  bxh = 0.5*bx


  if (nmp.shape(bx)[0] != frs):
    print "ERROR:  Number of frames in box variable does not match position variable\n"
    return 0



  # Pick log-spaced times to aid in fitting later on
  times = nmp.logspace(0,m.log10(frs), NP+3)
  delts = nmp.zeros(NP+3, 'i')

  for i in range(0, NP+3):
    delts[i] = int(times[i])
  delts = nmp.unique(delts)
  NP = nmp.shape(delts)[0]

  otp = open('to_c_dynamics.bin','wb')
  io.write_c_binary(otp,x,bx,delts)
  io.write_int(otp,q[0])
  io.write_int(otp,q[1])
  io.write_int(otp,q[2])
  otp.close()

  os.system('python_nonaff_dynamics to_c_dynamics.bin')

  # Read in output
  dat = nmp.zeros((NP,4),'d')
  inp = open('dynamics.out','r')
  for i in range(0,NP):
    line = inp.readline().split()
    dat[i,0] = float(line[0])
    dat[i,1] = float(line[1])
    dat[i,2] = float(line[2])
    dat[i,3] = float(line[3])
  inp.close()

  return dat

