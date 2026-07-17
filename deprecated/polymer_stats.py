import numpy as nmp
import math as m
import pbc_utils as pbc

def calc_Cinf(x, box, nch, chl, N):
  boxh = 0.5*box
  ns = nch*chl
  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i], box, boxh)

  cos_avg = 0.0
  ct = 0.0
  for i in range(0,nch):
    for j in range(0,chl-N-1):
      i1 = i*chl
      i2 = i1 + 1
      
      j1 = i1 + N
      j2 = j1 + 1

      dx12 = pbc.pbc_vdr(x[i1],x[i2],box,boxh)
      dx23 = pbc.pbc_vdr(x[j1],x[j2],box,boxh)

      dot = nmp.dot(dx12,dx23)
      mg1 = m.sqrt(nmp.dot(dx12,dx12))
      mg2 = m.sqrt(nmp.dot(dx23,dx23))

      cos_avg += (dot / (mg1 * mg2) )
      ct += 1.0

  cos_avg *= (1.0 / ct)

  Cinf =  (1+cos_avg) / (1-cos_avg)
  
  return Cinf

def connect_chains(x, box, nch, chl):
  boxh = 0.5*box
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = i*chl + j
      i2 = i1 + 1
      
      dr = pbc.pbc_vdr(x[i1],x[i2],box,boxh)
      
      x[i2] = x[i1] - dr

  return x


def avg_RNete(x, box, nch, chl, N):
  boxh = 0.5*box
  ns = nch*chl
  
  x = connect_chains(x, box, nch, chl)
  
  
  R2ete_avg = 0.0
  ct = 0.0
  for i in range(0,nch):
    for j in range(0, chl-N):
      i1 = i*chl + j
      i2 = i1 + N
    
      dx12 = x[i1] - x[i2]
      mag = nmp.dot(dx12,dx12)

      R2ete_avg += mag
      ct += 1.0

  return (R2ete_avg/ct)

def R2ete_all(x, box, nch, chl):
  boxh = 0.5*box
  ns = nch*chl
  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i], box, boxh)
  
  x = connect_chains(x, box, nch, chl)

  R2ete_list = []
  ct = 0.0
  for i in range(0,nch):
    i1 = i*chl
    i2 = i1 + chl-1
    
    #dx12 = pbc.pbc_vdr(x[i1],x[i2],box,boxh)
    dx12 = x[i1] - x[i2]
    mag = nmp.dot(dx12,dx12)

    R2ete_list.append(mag)

  return R2ete_list

def avg_R2ete(x, box, nch, chl):
  boxh = 0.5*box
  ns = nch*chl
  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i], box, boxh)
  
  x = connect_chains(x, box, nch, chl)

  R2ete_avg = 0.0
  ct = 0.0
  for i in range(0,nch):
    i1 = i*chl
    i2 = i1 + chl-1
    
    #dx12 = pbc.pbc_vdr(x[i1],x[i2],box,boxh)
    dx12 = x[i1] - x[i2]
    mag = nmp.dot(dx12,dx12)

    R2ete_avg += mag
    ct += 1.0

  return (R2ete_avg/ct)
