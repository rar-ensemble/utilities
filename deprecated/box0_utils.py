###########################################
# read_box0.py      6 May 2008
# Rob Riggleman
# Reads in a .box0 formatted file
# and returns the positions in 'x',
# labels in 'tp', box dimensions in 'box', 
# and number os particles in 'ns'
###########################################

import numpy as nmp

def read_all(name):
    inp = open(name,'r')

    L = inp.readline()
    L = inp.readline().split()
    strain = 0.0
    sh_drx = 0.0 
    if ( len( L ) > 1 ):
      strain = float( L[3] )
      if ( len( L ) > 4 ): 
        sh_drx = float( L[4] )

    L = inp.readline()
    L = L.split()
    ns = int(L[0])

    L = inp.readline()
    L = inp.readline()
    L = L.split()
    box = nmp.zeros(3,'d')
    box[0] = float(L[0])
    box[1] = float(L[1])
    box[2] = float(L[2])

    L = inp.readline()
    L = inp.readline()
    L = inp.readline()
    
    count = 0
 
    tp = nmp.zeros((ns,1),'c')
    x = nmp.zeros((ns,3),'d')

    for i in range(0,ns):

      L = inp.readline().split()
      if ( len(L) == 0):
        print "error on site",i

      tp[i] = (L[1])[0]
      x[i][0] = float(L[2])
      x[i][1] = float(L[3])
      x[i][2] = float(L[4])

      count += 1
      
    inp.close()

    if (strain != 0.0 ):
      print "Found a shear strain of", strain


      

    return x, tp, ns, box, strain, sh_drx

def write_all(ns, x, tp, box, name):
  inp = open(name,'w')

  inp.write('box0 file created by box0_utils.py\n\n')

  line = '%d sites in the file\n\n' % ns
  inp.write(line)

  line = '%f %f %f box dimensions\n\n' % (box[0],box[1],box[2])
  inp.write(line)

  inp.write('xyz positions of the particles\n\n')
  for i in range(0,ns):
    line = '%d %s %f %f %f\n' % (i,str(tp[i][0]),x[i][0],x[i][1],x[i][2])
    inp.write(line)

  inp.close()
