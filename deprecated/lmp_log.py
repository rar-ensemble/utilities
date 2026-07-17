import numpy as nmp

def read_log(name):
  inp = open(name,'r')

  L = ["stuff","stf"]
  lines = 0
  while ( L[0] != 'run'):
    L = inp.readline().split()
    if (len(L)==0):
      L = ['asdfasdf', 'jklewd']
    lines += 1

  L = inp.readline().split()
  L = inp.readline().split()
  print L
  
  ny = len(L)
  dat = nmp.zeros((10000, ny),'d')
  ndat = 0

  L = inp.readline().split()
  while (L[0] != 'Loop'):
    for j in range(0,ny):
      dat[ndat,j] = float(L[j])
    ndat += 1
    L = inp.readline().split()

  dat = dat[0:ndat,:]
  print 'shape of the data matrix:', nmp.shape(dat)

  inp.close()

  return dat

