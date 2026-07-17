import numpy as nmp

def write_frame( nm , x , tp ):
  ns = nmp.shape( x )[0]
  
  otp = open( nm , 'w' )
 
  line = '%d\n\n' % ns 
  otp.write(line )

  for i in range( 0 , ns ):
    line = '%s %f %f %f\n' % ( tp[i] , x[i,0] , x[i,1] , x[i,2] )
    otp.write(line)

  otp.close()

def read_frame(inp):
  # Read header
  L = inp.readline().split()
  
  if (len(L)==0):
    print "End of xyz file!\n"
    return 0

  ns = int(L[0])
  L = inp.readline()

  # Allocate positions
  x = nmp.zeros((ns,3),'double')
  tp = nmp.zeros(ns,'c')

  for i in range(0,ns):
    L = inp.readline().split()
    if (len(L)==4):
      tp[i] = L[0]
      x[i,0] = float(L[1])
      x[i,1] = float(L[2])
      x[i,2] = float(L[3])
    else:
      print "Error in xyz file!"
      return

  return [tp, x, ns]


# Read first frame only
def read_single(name):
  inp = open(name,'r')
  
  rt = read_frame(inp)

  inp.close()

  return rt


#  Reads in frame number "fr", counting from 1
def read_number(name, fr):
  inp = open(name,'r')

  for i in range(0,fr+1):
    rt = read_frame(inp)

  inp.close()
  print "Return contains [types, X, ns]"
  return rt


def read_trajectory(name):
  inp = open(name,'r')

  # Read the first frame
  rt = read_frame(inp)
  
  tp = rt[0]
  ns = rt[2]

  # Initially allocate 50 frames
  NFR = 50
  fr = 0
  x = nmp.zeros((NFR,ns,3),'d')

  
  while (rt!=0):
    x[fr,:,:] = rt[1]
    fr += 1

    if (fr==NFR):
      NFR += 50
      x = nmp.resize(x,(NFR,ns,3))
    
    rt = read_frame(inp)
  
  print "Read", fr, "frames in", name
  
  x = x[0:fr,:,:]
  inp.close()

  print "Return contains [tp, x, ns]\n"
  return [tp, x, ns]

