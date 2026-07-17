import numpy as nmp
import pbc_utils as pbc
import sys
import os
import random as rand


def write_uncharged( nm , index , mID , type , box , x , bonds , angles ):
  nangs = nmp.shape(angles)[0]
  nbonds = nmp.shape(bonds)[0]
  ns = nmp.shape(x)[0]

  otp = open( nm , 'w' )

  otp.write( 'arrrr\n\n' )

  ln = '%d atoms\n' % ns
  otp.write( ln )

  ln = '%d bonds\n' % nbonds
  otp.write( ln )

  ln = '%d angles\n\n' % nangs
  otp.write( ln )

  nang_types = 0
  if ( nangs > 0 ):
    nang_types = 1

  atom_types = nmp.shape(nmp.unique(type))[0]

  nbond_types = 0
  if ( nbonds > 0 ):
    nbond_types = 1

  ln = '%d atom types\n%d bond types\n%d angle types\n\n' % (atom_types,nbond_types,nang_types)
  otp.write( ln )

  ln = '%f %f xlo xhi\n' % (-box[0]/2.0 , box[0]/2.0)
  otp.write( ln )

  ln = '%f %f ylo yhi\n' % (-box[1]/2.0 , box[1]/2.0)
  otp.write( ln )

  ln = '%f %f zlo zhi\n\n' % (-box[2]/2.0 , box[2]/2.0)
  otp.write( ln )

  ln = 'Masses\n\n'
  otp.write( ln )
  for i in range(0, atom_types):
    ln = '%d 1.0\n' % (i+1)
    otp.write(ln)

  otp.write( '\nAtoms\n\n' )

  for i in range ( 0 , ns ):
    ln = '%d %d %d  %f %f %f\n' % ( index[i] , mID[i] , type[i] , x[i][0] , x[i][1] , x[i][2] )
    otp.write( ln )


  if ( nbonds > 0 ):
    otp.write( '\nBonds\n\n' )
    for i in range( 0 , nbonds ):
      ln = '%d %d  %d %d\n' % (i+1 , bonds[i][0] , bonds[i][1] , bonds[i][2])
      otp.write( ln )

  if ( nangs > 0 ):

    otp.write( 'Angles\n\n' )
    for i in range( 0 , nangs ):
      ln = '%d %d  %d %d %d\n' % ( i , angles[i][0] , angles[i][1] , angles[i][2] , angles[i][3])
      otp.write( ln )


  otp.close()

def write_charged( name , index, mID , type , charge , x , bx , bonds ):
  ns = nmp.shape( x ) [0]
  nbonds = nmp.shape( bonds ) [0]

  otp = open( name , 'w' )

  otp.write( 'huh?\n\n' )

  ln = '%d atoms\n' % ( ns )
  otp.write( ln )

  ln = '%d bonds\n\n' % ( nbonds )
  otp.write( ln )

  ln = '7 atom types\n1 bond types\n\n' 
  otp.write( ln )

  ln = '%f %f xlo xhi\n' % ( -bx[0]/2. , bx[0]/2. )
  otp.write(ln)

  ln = '%f %f ylo yhi\n' % ( -bx[1]/2. , bx[1]/2. )
  otp.write(ln)

  ln = '%f %f zlo zhi\n\n' % ( -bx[2]/2. , bx[2]/2. )
  otp.write(ln)


  otp.write('Masses\n\n')
  otp.write('1 1\n2 1\n3 1\n4 1\n5 1\n6 1\n7 1\n\n')

  otp.write( 'Atoms\n\n' )
  for i in range( 0 , ns ):
    ln = '%d %d %d  %1.2f  %f %f %f\n' % ( index[i] , mID[i], type[i], charge[i], x[i,0] , x[i,1] , x[i,2] )
    otp.write( ln )

  otp.write('\nBonds\n\n')
  for i in range( 0 , nbonds ):
    ln = '%d %d  %d %d\n' % ( i+1 , bonds[i,0] , bonds[i,1]+1 , bonds[i,2]+1 )
    otp.write( ln )

  otp.close()

def read_charged( name ):
  if (os.path.exists(name)==False):
    print "File", name, "not found!!"
    return 0.0

  ns = 1
  nbonds = 0
  nangles = 0
  box = nmp.zeros(3,'d')


  inp = open(name,'r')

  line = inp.readline().split()
  while (line != ""):
    if (len(line)==1 and line[0]=="Atoms"):
      break

    if (len(line)==2 and line[1] == "atoms"):
      ns = int(line[0])

    if (len(line)==2 and line[1] == "bonds"):
      nbonds = int(line[0])

    if (len(line)==2 and line[1] == "angles"):
      nangles = int(line[0])

    if (len(line)==4 and line[2] == "xlo"):
      box[0] = float(line[1]) - float(line[0])

    if (len(line)==4 and line[2] == "ylo"):
      box[1] = float(line[1]) - float(line[0])
    
    if (len(line)==4 and line[2] == "zlo"):
      box[2] = float(line[1]) - float(line[0])
    

    line = inp.readline().split()

  print "Read", ns, "atoms", nbonds, "bonds and", box, "box dimensions from header of", name

  #Allocate memory for bonds and atoms
  x = nmp.zeros((ns,3),'d')
  charge = nmp.zeros( ns , 'd' )
  type = nmp.zeros(ns,'i')
  mID = nmp.zeros(ns,'i')
  index = nmp.zeros(ns,'i')
  bonds = nmp.zeros((nbonds, 3),'i')
  angles = nmp.zeros((nangles,4),'i')

  line = inp.readline().split()

  for i in range(0,ns):
    line = inp.readline().split()

    if (len(line)<6):
      print "Error on line", i
      break

    id = int( line[0] ) - 1
    index[id] = int(line[0])
    mID[id] = int(line[1])
    type[id] = int(line[2])
    charge[id] = float( line[3] )
    x[id][0], x[id][1], x[id][2] = float(line[4]), float(line[5]), float(line[6])

  line = inp.readline().split()
  line = inp.readline().split()
  if (line==[]):
    print "bonds header not found!"
    print "Return contains [index, molecule-ID, type, x, box]\n"
    return [index, mID, type, charge, x, box]
  elif (line[0]=="Velocities"):
    for j in range(0,ns+3):
      inp.readline().split()

  line = inp.readline().split()
  for j in range(0,nbonds):
    line = inp.readline().split()
    if (len(line) < 4):
      print "Error on bond",j
      break

    bonds[j,0], bonds[j,1], bonds[j,2] = int(line[1]), int(line[2])-1, int(line[3])-1

  line = inp.readline().split()
  line = inp.readline().split()
  if (line==[]):
    print "Return contains [index, molecule-ID, type, x, box, bonds]\n"
    return [index, mID, type, charge , x, box, bonds ]

  line = inp.readline()
  for j in range(0,nangles):
    line = inp.readline().split()
    if (len(line) < 5):
      print "Error on angle",j
      break

    angles[j,0] = int(line[1])
    angles[j,1] = int(line[2])
    angles[j,2] = int(line[3])
    angles[j,3] = int(line[4])

  print "Return contains [index, molecule-ID, type, x, box, bonds, angles]\n"
  return [index, mID, type, charge, x, box, bonds, angles]



def get_is_energy(name):
  inp = open(name,'r')
  line = "asdfadsf"
  ct = 1
  while (inp):
    
    line = inp.readline().split()
    print line
    if (len(line) > 0 and line[0]=="Energy"):
      break
    
    ct += 1
    if (ct==100):
      break

  line = inp.readline().split()
  inp.close()

  return float(line[2])

    
def write_lammps_swap(x, box, npar, nch, chl, name):
  ns = nmp.shape(x)[0]

  if (ns != (npar*50 + nch*chl)):
    print 'Error in the number of sites!'
    return

  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (nch*chl + npar*50)
  otp.write(line)

  line = '%d bonds\n' % (nch*(chl-1) + 49*npar)
  otp.write(line)
  otp.write('0 angles\n')
  otp.write('0 dihedrals\n')
  otp.write('0 impropers\n\n')

  if (npar > 0):
    otp.write('3 atom types\n2 bond types\n')
  else:
    otp.write('1 atom types\n1 bond types\n')
  otp.write('0 angle types\n')
  otp.write('0 dihedral types\n')
  otp.write('0 improper types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  if (npar > 0):
    otp.write('2 1.000000\n')
    otp.write('3 1.000000\n')
  otp.write('Atoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      mid = -1 ;
      if ( j < chl / 2 ):
        mid = j + 1
      else:
        mid = chl - j
      line = '%d %d %d     %f  %f  %f\n' % (i1,mid,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  # Write the particle positions
  for i in range(nch, nch+npar):
    line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,3, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1
    for j in range(1,50):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1
  # Write the particle bonds
  for i in range(0,npar):
    i1 = nch*chl + i*50 + 1
    for j in range(1,50):
      i2 = i1 + j
      line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
      otp.write(line)
      bi += 1

  otp.close()


    
def write_lammps_mg_surf(x, nsurf, box, nmol,  name):
  ns = nmp.shape(x)[0]


  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (ns)
  otp.write(line)

  line = '%d bonds\n' % ( 6 * nmol )
  otp.write(line)

  line = '%d angles\n' % ( 6 * nmol )
  otp.write(line)


  otp.write('2 atom types\n2 bond types\n')
  otp.write('2 angle types\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)

  line = '%f   %f zlo zhi\n' % ( nmp.min(x[:,2])-1.0 , nmp.max(x[:,2])+5.0 )
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  otp.write('2  1.000000\n')
  otp.write('\nAtoms\n\n')
  
  
  i1 = 1
  for i in range( 0 , nsurf ):
    line = '%d %d %d     %f  %f  %f\n' % (i1,1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1

  for i in range(0,nmol):
    for j in range(0,7):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+2,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nmol):
    io = 7 * i + nsurf 
    i1 = io + 1
    i2 = io + 2
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 3
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 4
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 4
    i2 = io + 5
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 6
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 6
    i2 = io + 7
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1


  line = '\nAngles\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer angle
  for i in range(0,nmol):
    io = 7 * i + nsurf 
    i1 = io + 1
    i2 = io + 2
    i3 = io + 3
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 4
    i3 = io + 5
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 6
    i3 = io + 7
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 1
    i3 = io + 4
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 1
    i3 = io + 6
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 4
    i2 = io + 1
    i3 = io + 6
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1



  otp.close()

    
def write_lammps_mg_tamd(x, box, nmol,  name):
  ns = nmp.shape(x)[0]


  otp = open(name,'w')


  otp.write('Generated by RAR python library\n\n')

  line = '%d atoms\n' % (ns)
  otp.write(line)

  line = '%d bonds\n' % ( 6 * nmol + nmol*7)
  otp.write(line)

  line = '%d angles\n' % ( 6 * nmol )
  otp.write(line)

  otp.write('0 dihedrals\n')
  otp.write('0 impropers\n\n')

  otp.write('2 atom types\n3 bond types\n')
  otp.write('2 angle types\n')
  otp.write('0 dihedral types\n')
  otp.write('0 improper types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  otp.write('2  1.000000\n')
  otp.write('\nAtoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  for i in range(0,nmol):
    for j in range(0,7):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  for i in range(0,nmol):
    for j in range(0,7):
      line = '%d %d %d     %f  %f  %f\n' % (i1,nmol+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nmol):
    io = 7 * i
    i1 = io + 1
    i2 = io + 2
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 3
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 4
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 4
    i2 = io + 5
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 6
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 6
    i2 = io + 7
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

  for i in range(0, nmol*7):
    i1 = i+1
    i2 = i1 + nmol*7
    line = '%d 3  %d %d\n' %(bi, i1, i2)
    otp.write(line)
    bi += 1

  line = '\nAngles\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nmol):
    io = 7 * i
    i1 = io + 1
    i2 = io + 2
    i3 = io + 3
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 4
    i3 = io + 5
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 6
    i3 = io + 7
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 1
    i3 = io + 4
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 1
    i3 = io + 6
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 4
    i2 = io + 1
    i3 = io + 6
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1



  otp.close()
    
def write_lammps_mg(x, box, nmol,  name):
  ns = nmp.shape(x)[0]


  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (ns)
  otp.write(line)

  line = '%d bonds\n' % ( 6 * nmol )
  otp.write(line)

  line = '%d angles\n' % ( 6 * nmol )
  otp.write(line)

  otp.write('0 dihedrals\n')
  otp.write('0 impropers\n\n')

  otp.write('1 atom types\n2 bond types\n')
  otp.write('2 angle types\n')
  otp.write('0 dihedral types\n')
  otp.write('0 improper types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  otp.write('\nAtoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  for i in range(0,nmol):
    for j in range(0,7):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nmol):
    io = 7 * i
    i1 = io + 1
    i2 = io + 2
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 3
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 4
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 4
    i2 = io + 5
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 6
    line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
    otp.write(line)
    bi += 1

    i1 = io + 6
    i2 = io + 7
    line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
    otp.write(line)
    bi += 1


  line = '\nAngles\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nmol):
    io = 7 * i
    i1 = io + 1
    i2 = io + 2
    i3 = io + 3
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 4
    i3 = io + 5
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 1
    i2 = io + 6
    i3 = io + 7
    line = '%d %d  %d %d %d\n' % (bi, 1, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 1
    i3 = io + 4
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 2
    i2 = io + 1
    i3 = io + 6
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1

    i1 = io + 4
    i2 = io + 1
    i3 = io + 6
    line = '%d %d  %d %d %d\n' % (bi, 2, i1, i2, i3)
    otp.write(line)
    bi += 1



  otp.close()

    
def write_lammps_supported(x, box, ns_slab, nch, chl, name):
  ns = nmp.shape(x)[0]

  if (ns != (ns_slab + nch*chl)):
    print 'Error in the number of sites!'
    print 'ns =', ns, 'ns_slab + nch*chl =', ns_slab + nch*chl
    return

  otp = open(name,'w')


  otp.write('Generated by RARs messy code\n\n')

  line = '%d atoms\n' % (ns)
  otp.write(line)

  line = '%d bonds\n' % (nch*(chl-1))
  otp.write(line)

  if (ns_slab > 0):
    otp.write('\n2 atom types\n1 bond types\n\n')
  else:
    otp.write('\n1 atom types\n1 bond types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1 1.000000\n')
  if (ns_slab > 0):
    otp.write('2 1.000000\n')
  otp.write('\nAtoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the slab positions
  for i in range(0,ns_slab):
    line = '%d %d %d     %f  %f  %f\n' % (i1,1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1
  
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+2,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1 + ns_slab
      i2 = i1 + 1 
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1

  otp.close()

    
def write_lammps_2d_mix(x, box, fracA, nch, chl, name):
  ns = nmp.shape(x)[0]
  npar = 0

  if (ns != (npar*50 + nch*chl)):
    print 'Error in the number of sites!'
    return

  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (nch*chl + npar*50)
  otp.write(line)

  line = '%d bonds\n' % (nch*(chl-1) + 49*npar)
  otp.write(line)

  otp.write('2 atom types\n0 bond types\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  otp.write('2 1.000000\n')
  otp.write('Atoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      if ( rand.random() < fracA ):
        line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      else:
        line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  # Write the particle positions
  for i in range(nch, nch+npar):
    line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,3, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1
    for j in range(1,50):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1
  # Write the particle bonds
  for i in range(0,npar):
    i1 = nch*chl + i*50 + 1
    for j in range(1,50):
      i2 = i1 + j
      line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
      otp.write(line)
      bi += 1

  otp.close()


    
def write_lammps_angles(x, box, nsolv, nch, chl, name):
  ns = nmp.shape(x)[0]

  if (ns != (nsolv + nch*chl)):
    print 'Error in the number of sites!'
    return

  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (nch*chl + nsolv)
  otp.write(line)

  line = '%d bonds\n' % (nch*(chl-1) )
  otp.write(line)
  line = '%d angles\n\n' % (nch*(chl-2))
  otp.write(line)

  otp.write('2 atom types\n1 bond types\n')
  otp.write('1 angle types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  otp.write('2  1.000000\n')
  otp.write('Atoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  # Write the particle positions
  for i in range(0, nsolv):
    line = '%d %d %d     %f  %f  %f\n' % (i1,nch+i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1

  line = '\nAngles\n\n'
  otp.write(line)
  ai = 1
  # Write the polymer angles
  for i in range(0,nch):
    for j in range(0,chl-2):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      i3 = i1 + 2
      line = '%d %d   %d %d %d\n' % (ai, 1, i1, i2, i3)
      otp.write(line)
      ai += 1


  otp.close()


    
def write_lammps(x, box, npar, nch, chl, name):
  ns = nmp.shape(x)[0]

  if (ns != (npar*50 + nch*chl)):
    print 'Error in the number of sites!'
    return

  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (nch*chl + npar*50)
  otp.write(line)

  line = '%d bonds\n' % (nch*(chl-1) + 49*npar)
  otp.write(line)
  otp.write('0 angles\n')
  otp.write('0 dihedrals\n')
  otp.write('0 impropers\n\n')

  if (npar > 0):
    otp.write('3 atom types\n2 bond types\n')
  else:
    otp.write('1 atom types\n1 bond types\n')
  otp.write('0 angle types\n')
  otp.write('0 dihedral types\n')
  otp.write('0 improper types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  if (npar > 0):
    otp.write('2 1.000000\n')
    otp.write('3 1.000000\n')
  otp.write('Atoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  # Write the particle positions
  for i in range(nch, nch+npar):
    line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,3, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1
    for j in range(1,50):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1
  # Write the particle bonds
  for i in range(0,npar):
    i1 = nch*chl + i*50 + 1
    for j in range(1,50):
      i2 = i1 + j
      line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
      otp.write(line)
      bi += 1

  otp.close()



def write_antiplas(x, box, npar, nch, chl, name):
  ns = nmp.shape(x)[0]

  otp = open(name,'w')


  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (nch*chl + npar)
  otp.write(line)

  line = '%d bonds\n' % (nch*(chl-1) )
  otp.write(line)
  otp.write('0 angles\n')
  otp.write('0 dihedrals\n')
  otp.write('0 impropers\n\n')

  if (npar > 0):
    otp.write('2 atom types\n1 bond types\n')
  else:
    otp.write('1 atom types\n1 bond types\n')
  otp.write('0 angle types\n')
  otp.write('0 dihedral types\n')
  otp.write('0 improper types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  if (npar > 0):
    otp.write('2 0.125000\n')
  otp.write('Atoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      mid = 0
      if (j < chl/2):
        mid = j+1
      elif (j==(chl/2)):
        mid = chl/2
      else:
        mid = chl - j
      #line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      line = '%d %d %d     %f  %f  %f\n' % (i1,mid,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  # Write the particle positions
  for i in range(nch, nch+npar):
    line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
    otp.write(line)
    i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1

  otp.close()

def read_log_dat(name):
  inp = open(name,'r')

  save_fr = 0
  total_steps = 0
  np = 2

  L = inp.readline().split()
  while (inp):
    L = inp.readline().split()
    if (len(L) > 1 and L[0] == 'thermo'):
      save_fr = int(L[1])

    if (len(L) > 1 and L[0] == 'run'):
      total_steps = int(L[1])

    if (len(L) > 0 and L[0] == 'Step'):
      break

  ny = len(L)
  if (save_fr>0):
    np = total_steps/save_fr + 1

  print "Reading in the following prop [",ny, ',', np, ']\n', L

  dat = nmp.zeros((np, ny),'d')
  for i in range(0,np):
    L = inp.readline().split()
    for j in range(0, len(L) ):
      dat[i,j] = float(L[j])

  inp.close()

  return dat




def read_data(name):
  if (os.path.exists(name)==False):
    print "File", name, "not found!!"
    return 0.0

  ns = 1
  nbonds = 0
  nangles = 0
  box = nmp.zeros(3,'d')


  inp = open(name,'r')

  line = inp.readline().split()
  while (line != ""):
    if (len(line) >= 1 and line[0]=="Atoms"):
      break

    if (len(line)==2 and line[1] == "atoms"):
      ns = int(line[0])

    if (len(line)==2 and line[1] == "bonds"):
      nbonds = int(line[0])

    if (len(line)==2 and line[1] == "angles"):
      nangles = int(line[0])

    if (len(line)==4 and line[2] == "xlo"):
      box[0] = float(line[1]) - float(line[0])

    if (len(line)==4 and line[2] == "ylo"):
      box[1] = float(line[1]) - float(line[0])
    
    if (len(line)==4 and line[2] == "zlo"):
      box[2] = float(line[1]) - float(line[0])
    

    line = inp.readline().split()

  print "Read", ns, "atoms", nbonds, "bonds and", box, "box dimensions from header of", name

  #Allocate memory for bonds and atoms
  x = nmp.zeros((ns,3),'d')
  type = nmp.zeros(ns,'i')
  mID = nmp.zeros(ns,'i')
  index = nmp.zeros(ns,'i')
  bonds = nmp.zeros((nbonds, 3),'i')
  angles = nmp.zeros((nangles,4),'i')

  line = inp.readline().split()

  for i in range(0,ns):
    line = inp.readline().split()

    if (len(line)<6):
      print "Error on line", i
      break

    ind = int( line[0] ) - 1
    index[ind] = ind+1
    mID[ind] = int(line[1])
    type[ind] = int(line[2])
    x[ind][0], x[ind][1], x[ind][2] = float(line[3]), float(line[4]), float(line[5])

  line = inp.readline().split()
  line = inp.readline().split()
  if (line==[]):
    print "bonds header not found!"
    print "Return contains [index, molecule-ID, type, x, box]\n"
    return [index, mID, type, x, box]
  elif (line[0]=="Velocities"):
    for j in range(0,ns+3):
      inp.readline().split()

  line = inp.readline().split()
  for j in range(0,nbonds):
    line = inp.readline().split()
    if (len(line) < 4):
      print "Error on bond",j
      break

    bonds[j,0], bonds[j,1], bonds[j,2] = int(line[1]), int(line[2]), int(line[3])

  line = inp.readline().split()
  line = inp.readline().split()
  if (line==[]):
    print "Return contains [index, molecule-ID, type, x, box, bonds]\n"
    return [index, mID, type, x, box, bonds]

  line = inp.readline()
  for j in range(0,nangles):
    line = inp.readline().split()
    if (len(line) < 5):
      print "Error on angle",j
      break

    angles[j,0] = int(line[1])
    angles[j,1] = int(line[2])
    angles[j,2] = int(line[3])
    angles[j,3] = int(line[4])

  print "Return contains [index, molecule-ID, type, x, box, bonds, angles]\n"
  return [index, mID, type, x, box, bonds, angles]


def deshuffle_bonds(bonds):
  nbonds = nmp.shape(bonds)[0]
  ind = 1
  rebond = nmp.zeros((nbonds,3),'i')
  old_bond_flag = nmp.zeros(nbonds,'i')
  ns = nmp.max(bonds[:,1:3])
  ends = nmp.zeros(ns+1,'i')
  particle_bonds = nmp.zeros(ns+1,'i')
  
  # Number of bonds on each site
  for j in range(0,nbonds):
    id1 = bonds[j,1]
    id2 = bonds[j,2]
    particle_bonds[id1] += 1
    particle_bonds[id2] += 1

 
  list_3 = nmp.where( particle_bonds[:] > 2 )
  n3 = len( list_3 )  

  if ( n3 > 1 or list_3[0] != [] ):
    print "Something screwy happened", n3, "particles have 3 bonds!"
    for i in range( 0 , n3 ):
      print list_3[i]

  # Chain/particle ends
  ed = nmp.where(particle_bonds==1)

  print 'Chain ends\n',ed

  # Flag marking chain/particle ends
  ends[ed] = 1

  # Find first bond and go from there
  nnew = 0
  for j in range(0,nbonds):
    if (bonds[j,1]==1):
      old_bond_flag[j] = 1
      rebond[0,0:3] = bonds[j,:]
      nnew += 1

  tries = 0
  while (nnew < nbonds):
    id1 = rebond[nnew-1,1]
    id2 = rebond[nnew-1,2]

    found_flag=0
    for j in range(0,nbonds):
      if ( (bonds[j,1]==id1 and bonds[j,2]!=id2) or (bonds[j,1]==id2 and bonds[j,2]!=id1) ):
        rebond[nnew,0:3] = bonds[j,0:3]
        old_bond_flag[j] = 1
        found_flag=1
        nnew += 1
        break

    if (found_flag==0):
      for j in range(0,nbonds):
        id1 = bonds[j,1]
        id2 = bonds[j,2]
        if (old_bond_flag[j]==0 and (ends[id1]==1 or ends[id2==1]) ):
          old_bond_flag[j] = 1
          rebond[nnew,0:3] = bonds[j,0:3]
          rebond[nnew,2] = bonds[j,2]
          nnew += 1
          break


    tries += 1
    if (tries == 5*nbonds):
      print "error, giving up!\n"

    if (nnew % 500 == 0):
      print "Found" , nnew , "of", nbonds

  si= nmp.zeros(ns,'i')
  si[0] = rebond[0,1]
  si[1] = rebond[0,2]
  ind = 2
  for j in range(1,nbonds):
    if (si[ind-1]==rebond[j,1]):
      si[ind] = rebond[j,2]
      ind += 1
    else:
      si[ind] = rebond[j,1]
      ind += 1
      si[ind] = rebond[j,2]
      ind += 1

  print ind, ns
  print "Return contains [sorted bonds, site index]"
  return [rebond, si]



def write_smooth(x, box, name):
  ns = nmp.shape(x)[0]

  nch = 46
  chl = 500


  print "Enter number of sites per nanorod\n"

  inp = sys.__stdin__
  line = inp.readline().split()

  rodL = int(line[0])
  nrod = (ns - nch*chl) / rodL

  otp = open(name,'w')

  otp.write('Generated by RARs stupid code\n\n')

  line = '%d atoms\n' % (ns)
  otp.write(line)

  line = '%d bonds\n' % (nrod*(rodL-1) + nch*(chl-1) )
  otp.write(line)
  line = '%d angles\n' % (nrod*(rodL-2))
  if (rodL==1):
    line = '0 angles\n'
  otp.write(line)

  otp.write('0 dihedrals\n')
  otp.write('0 impropers\n\n')

  if (nrod > 0):
    otp.write('2 atom types\n2 bond types\n1 angle types\n')
  else:
    otp.write('1 atom types\n1 bond types\n0 angle types\n')
  otp.write('0 dihedral types\n')
  otp.write('0 improper types\n\n')

  boxh = 0.5*box
  line = '%f   %f xlo xhi\n' % (-boxh[0], boxh[0])
  otp.write(line)
  line = '%f   %f ylo yhi\n' % (-boxh[1], boxh[1])
  otp.write(line)
  line = '%f   %f zlo zhi\n' % (-boxh[2], boxh[2])
  otp.write(line)

  otp.write('\nMasses\n\n')
  otp.write('1  1.000000\n')
  if (nrod > 0):
    otp.write('2 1.000000\n')
  otp.write('Atoms\n\n')
  

  for i in range(0,ns):
    x[i] = pbc.pbc_inbox(x[i],box,boxh)
  
  i1 = 1
  # Write the polymer positions
  for i in range(0,nch):
    for j in range(0,chl):
      mid = 0
      if (j < chl/2):
        mid = j+1
      elif (j==(chl/2)):
        mid = chl/2
      else:
        mid = chl - j
      #line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      line = '%d %d %d     %f  %f  %f\n' % (i1,mid,1, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  # Write the particle positions
  for i in range(nch, nch+nrod):
    for j in range(0,rodL):
      line = '%d %d %d     %f  %f  %f\n' % (i1,i+1,2, x[i1-1,0], x[i1-1,1], x[i1-1,2])
      otp.write(line)
      i1 += 1

  line = '\nBonds\n\n'
  otp.write(line)
  bi = 1

  # Write the polymer bonds
  for i in range(0,nch):
    for j in range(0,chl-1):
      i1 = chl*i + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 1, i1, i2)
      otp.write(line)
      bi += 1
  # Write the particle bonds
  for i in range(0,nrod):
    for j in range(0,rodL-1):
      i1 = nch*chl + i*rodL + j + 1
      i2 = i1 + 1
      line = '%d  %d  %d  %d\n' % (bi, 2, i1, i2)
      otp.write(line)
      bi += 1

  ai = 1
  
  if (rodL > 1):
    line = '\nAngles\n\n'
    otp.write(line)
    for i in range(0, nrod):
      for j in range(0, rodL-2):
        i1 = nch*chl + i*rodL + j + 1
        i2 = i1 + 1
        i3 = i2 + 1
 
        line = '%d  %d  %d  %d  %d\n' % (ai, 1, i1, i2, i3)
        otp.write(line)
        ai += 1

  otp.close()

def write_arl( nm , index , mID , type , box , x , bonds , angles ):
  nangs = nmp.shape(angles)[0]
  if ( angles == [] or angles == [[]] ):
    print "No angles, you're the angles!"
    nangs = 0

  nbonds = nmp.shape(bonds)[0]
  ns = nmp.shape(x)[0]

  otp = open( nm , 'w' )

  otp.write( 'arrrr\n\n' )

  ln = '%d atoms\n' % ns
  otp.write( ln )

  ln = '%d bonds\n' % nbonds
  otp.write( ln )

  ln = '%d angles\n\n' % nangs
  otp.write( ln )

  ln = '%d atom types\n%d bond types\n' % (nmp.shape( nmp.unique(type))[0], \
          nmp.shape(nmp.unique(bonds[:,0]))[0])
  otp.write( ln )
  
  print nmp.unique(angles[:,0])

  if ( nangs > 0 ):
    ln = '%d angle types\n\n' % ( nmp.shape(nmp.unique(angles[:,0]))[0] )
  else:
    ln = '0 angle types\n\n' 
  otp.write( ln )

  ln = '%f %f xlo xhi\n' % (-box[0]/2.0 , box[0]/2.0)
  otp.write( ln )

  ln = '%f %f ylo yhi\n' % (-box[1]/2.0 , box[1]/2.0)
  otp.write( ln )

  ln = '%f %f zlo zhi\n\n' % (-box[2]/2.0 , box[2]/2.0)
  otp.write( ln )

  ln = 'Masses\n\n1 1\n2 1\n\n'
  otp.write( ln )

  otp.write( 'Atoms\n\n' )

  for i in range ( 0 , ns ):
    ln = '%d %d %d  %f %f %f\n' % ( index[i] , mID[i] , type[i] , x[i][0] , x[i][1] , x[i][2] )
    otp.write( ln )


  otp.write( '\nBonds\n\n' )
  for i in range( 0 , nbonds ):
    ln = '%d %d  %d %d\n' % (i+1 , bonds[i][0] , bonds[i][1] , bonds[i][2])
    otp.write( ln )


  otp.write( '\nAngles\n\n' )
  for i in range( 0 , nangs ):
    ln = '%d %d  %d %d %d\n' % ( i+1 , angles[i][0] , angles[i][1] , angles[i][2] , angles[i][3])
    otp.write( ln )


  otp.close()




def write_cuzr( nm , index , type , box , x  ):
  ns = nmp.shape(x)[0]

  otp = open( nm , 'w' )

  otp.write( 'arrrr\n\n' )

  ln = '%d atoms\n' % ns
  otp.write( ln )



  ln = '%d atom types\n\n' % (2)
  otp.write( ln )

  ln = '%f %f xlo xhi\n' % (-box[0]/2.0 , box[0]/2.0)
  otp.write( ln )

  ln = '%f %f ylo yhi\n' % (-box[1]/2.0 , box[1]/2.0)
  otp.write( ln )

  ln = '%f %f zlo zhi\n\n' % (-box[2]/2.0 , box[2]/2.0)
  otp.write( ln )

  ln = 'Masses\n\n1 63.546\n2 91.224\n\n'
  otp.write( ln )

  otp.write( 'Atoms\n\n' )

  for i in range ( 0 , ns ):
    ln = '%d %d  %f %f %f\n' % ( index[i]  , type[i] , x[i][0] , x[i][1] , x[i][2] )
    otp.write( ln )


  otp.close()





















