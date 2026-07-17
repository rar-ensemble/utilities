import struct
import numpy as nmp

def append_file(nm , lst):
  otp = open(nm, 'a' )
  for i in range( 0 , len(lst) ):
    otp.write( lst[i] )
  otp.close()

def write_float(otp, D):
  d = struct.pack('d',D)
  otp.write(d)

def write_int(otp, N):
  d = struct.pack('i',N)
  otp.write(d)

def write_float_matrix(otp, D, nx, ny):
  for i in range(0,nx):
    for j in range(0,ny):
      write_float(otp,D[i][j])

def write_int_matrix(otp, N, nx, ny):
  for i in range(0,nx):
    for j in range(0,ny):
      write_int(otp,N[i][j])

def write_int_vector(otp, N, nx):
  for i in range(0,nx):
    write_int(otp,N[i])

def write_float_vector(otp, N, nx):
  for i in range(0,nx):
    write_float(otp,N[i])

def read_float(inp):
  tp = inp.read(8)
  D = struct.unpack('d',tp)
  D = float(D[0])
  return D

def read_int(inp):
  tp = inp.read(4)
  N = struct.unpack('i',tp)
  N = int(N[0])
  return N

def read_int_vector(inp,nx):
  A = nmp.zeros(nx,'i')
  for i in range(0,nx):
    A[i] = read_int(inp)
  return A

def read_int_matrix(inp,nx,ny):
  A = nmp.zeros((nx,ny),'i')

  for i in range(0,nx):
    for j in range(0,ny):
      A[i][j] = read_int(inp)
  return A

def read_float_matrix(inp,nx,ny):
  A = nmp.zeros((nx,ny),'d')

  for i in range(0,nx):
    for j in range(0,ny):
      A[i][j] = read_float(inp)
  return A

def read_vec_text(name, nx, ny):
  inp = open(name,'r')
  dat = nmp.zeros((ny,nx),'d')
  for i in range(0,ny):
    L = inp.readline().split()
    dat[i,0] = float(L[0])
    dat[i,1] = float(L[1])
  inp.close()
  return dat

def write_vec_text(dat,name):
  otp = open(name,'w')
  ny = nmp.shape(dat)[0]
  nx = nmp.shape(dat)[1]
  if (nx != 2):
    otp.write('dat has wrong shape!\n')
    otp.close()
    return
  for i in range(0,ny):
    line = '%e %e\n' % (dat[i][0],dat[i][1])
    otp.write(line)
  otp.close()
    

def save_state(MAX,Nstates,ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd):
  otp = open('saved-states.bin','wb')
  
  write_int(otp,MAX)
  write_int(otp,Nstates)
  write_int(otp,ns)
  write_float_matrix(otp,DATA,MAX,5)
  for i in range(0,MAX):
    write_float_matrix(otp,STATES[i],ns,3)
  for i in range(0,MAX):
    write_float_matrix(otp,SADDLE[i],ns,3)
  write_float_matrix(otp,SADDAT,MAX,4)
  write_float_matrix(otp,Mcon,MAX,MAX)
  write_int(otp,NE)
  write_int_vector(otp,ActInd,MAX)

  otp.close()

def read_state():
  inp = open('saved-states.bin','rb')
  MAX = read_int(inp)
  Nstates = read_int(inp)
  ns = read_int(inp)
  DATA = read_float_matrix(inp,MAX,5)
  
  STATES = nmp.zeros((MAX,ns,3),'d')
  for i in range(0,MAX):
    STATES[i] = read_float_matrix(inp,ns,3)
  
  SADDLE = nmp.zeros((MAX,ns,3),'d')
  for i in range(0,MAX):
    SADDLE[i] = read_float_matrix(inp,ns,3)

  SADDAT = read_float_matrix(inp,MAX,4)
  Mcon = read_float_matrix(inp,MAX,MAX)
  NE = read_int(inp)
  ActInd = read_int_vector(inp,MAX)
  inp.close()

  return MAX, Nstates, ns, DATA, STATES, SADDLE, SADDAT, Mcon, NE, ActInd


def resume():
  inp = open('resume.bin','r')
  MAX = read_int(inp)
  Nstates = read_int(inp)
  ns = read_int(inp)
  DATA = read_float_matrix(inp,MAX,5)
  
  STATES = nmp.zeros((MAX,ns,3),'d')
  for i in range(0,MAX):
    STATES[i] = read_float_matrix(inp,ns,3)
  
  SADDLE = nmp.zeros((MAX,ns,3),'d')
  for i in range(0,MAX):
    SADDLE[i] = read_float_matrix(inp,ns,3)

  SADDAT = read_float_matrix(inp,MAX,4)
  Mcon = read_float_matrix(inp,MAX,MAX)
  NE = read_int(inp)
  ActInd = read_int_vector(inp,MAX)
  inp.close()

  return MAX, STATES, SADDLE, SADDAT, DATA, Mcon, NE, Nstates, ActInd

def write_c_binary(otp, x, bx, delts):
  frs = nmp.shape(x)[0]
  ns = nmp.shape(x)[1]
  NP = nmp.shape(delts)[0]
  write_int(otp,frs)
  write_int(otp,ns)
  write_int(otp,NP)
  for t in range(0,frs):
        write_float_matrix(otp,x[t,:,:],ns,3)
  write_float_matrix(otp,bx,frs,3)
  write_int_vector(otp,delts,NP)

