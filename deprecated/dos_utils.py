import math as m
import cmath as cm
import numpy as nmp
def write_modA_3D_AVG( Om , H ):
  nx = nmp.shape(H)[0]
  ny = nmp.shape(H)[1]
  nz = nmp.shape(H)[2]
  otp = open( "avg_combined.dat" , "w" )

  for i in range( 0 , nx ):
    for j in range( 0 , ny ):
      for k in range( 0 , nz ):
        line = "%f %f %f %1.14e %1.14e\n" % (H[i,j,k,0] , H[i,j,k,1] , H[i,j,k,2] , Om[i,j,k].real , Om[i,j,k].imag)
        otp.write(line)

  otp.close()

def write_modA_3D_DOS( Om , H ):
  nx = nmp.shape(H)[0]
  ny = nmp.shape(H)[1]
  nz = nmp.shape(H)[2]
  otp = open( "dos_combined.dat" , "w" )

  for i in range( 0 , nx ):
    for j in range( 0 , ny ):
      for k in range( 0 , nz ):
        line = "%f %f %f %1.14e\n" % (H[i,j,k,0] , H[i,j,k,1] , H[i,j,k,2] , Om[i,j,k])
        otp.write(line)

  otp.close()


def combine_modA_two_AVG(Om1 , Om2 , H1 , H2 ):
  nx1 = nmp.shape(H1)[0]
  nx2 = nmp.shape(H2)[0]
  ny = nmp.shape(H1)[1]
  
  nz = nmp.shape(Om1)[2]


  IRmin = nmp.min(H1[:,:,:,0])
  IRmax = nmp.max(H2[:,:,:,0])
  dIR = H1[1,0,0,0] - H1[0,0,0,0]

  nIRnew = int( (IRmax-IRmin)/dIR ) + 1

  Omnew = nmp.zeros([nIRnew,ny,nz],'cdouble')
  Hnew = nmp.zeros([nIRnew, ny,nz, 3],'d')
  for i in range(0,nIRnew):
    Hnew[i,:,:,0] = IRmin + dIR*i
    Hnew[i,:,:,1] = H1[nx1-1,:,:,1]
    Hnew[i,:,:,2] = H1[nx1-1,:,:,2]

  total_n_common = 0

  for i in range(0,nIRnew):
    
    ncom = 0
    for j in range(0,nx1):
      if (Hnew[i,0,0,0] == H1[j,0,0,0]):
        Omnew[i,:,:] += Om1[j,:,:]
        ncom += 1
        break

    for j in range(0,nx2):
      if (Hnew[i,0,0,0] == H2[j,0,0,0]):
        Omnew[i,:,:] += Om2[j,:,:]
        ncom += 1
        break

    if (ncom > 1):
      Omnew[i,:,:] *= 0.5
      total_n_common += 1
  
  print total_n_common , 'points in common'

  return Omnew, Hnew




def combine_modA_two_DOS(Om1 , Om2 , H1 , H2 ):
  nx1 = nmp.shape(H1)[0]
  nx2 = nmp.shape(H2)[0]
  ny = nmp.shape(H1)[1]
  
  nz = nmp.shape(Om1)[2]

  Oshift = 0.0
  ncommon = 0
  # Find the average difference between commonly sampled H-space
  for i in range(0,nx1):
    for j in range(0,ny):
      
      for k in range(0,nx2):
        for m in range(0,ny):
          for n in range(0,nz):
            
            if (H1[i,j,n,0]==H2[k,m,n,0] and H1[i,j,n,1]==H2[k,m,n,1]
                and H1[i,j,n,2] == H2[k,m,n,2] ):
              
              if (Om2[k,m,n] > 0 and Om1[i,j,n] > 0):
                Oshift += Om2[k,m,n] - Om1[i,j,n]
                ncommon += 1.0

  if ( ncommon == 0 ):
    print "No points in common!"
    exit(1)

  Oshift *= (1.0 / ncommon)
  print "Oshift:" , Oshift
  
  # Perform the shifting
  for i in range(0,nx2):
    for j in range(0,ny):
      for k in range(0,nz):
        if (Om2[i,j,k] > 0.0):
          Om2[i,j,k] -= Oshift

  IRmin = nmp.min(H1[:,:,:,0])
  IRmax = nmp.max(H2[:,:,:,0])
  dIR = H1[1,0,0,0] - H1[0,0,0,0]

  nIRnew = int( (IRmax-IRmin)/dIR ) + 1

  Omnew = nmp.zeros([nIRnew,ny,nz],'d')
  Hnew = nmp.zeros([nIRnew, ny,nz, 3],'d')
  for i in range(0,nIRnew):
    Hnew[i,:,:,0] = IRmin + dIR*i
    Hnew[i,:,:,1] = H1[nx1-1,:,:,1]
    Hnew[i,:,:,2] = H1[nx1-1,:,:,2]

  for i in range(0,nIRnew):
    
    ncom = 0
    for j in range(0,nx1):
      if (Hnew[i,0,0,0] == H1[j,0,0,0]):
        Omnew[i,:,:] += Om1[j,:,:]
        ncom += 1
        break

    for j in range(0,nx2):
      if (Hnew[i,0,0,0] == H2[j,0,0,0]):
        Omnew[i,:,:] += Om2[j,:,:]
        ncom += 1
        break

    if (ncom > 1):
      Omnew[i,:,:] *= 0.5

  return Omnew, Hnew


def combine_two_DOS(Om1, Om2, H1, H2):
  nx1 = nmp.shape(H1)[0]
  nx2 = nmp.shape(H2)[0]
  ny = nmp.shape(H1)[1]
  
  nz = nmp.shape(Om1)[2]

  Oshift = 0.0
  ncommon = 0
  # Find the average difference between commonly sampled H-space
  for i in range(0,nx1):
    for j in range(0,ny):
      
      for k in range(0,nx2):
        for m in range(0,ny):
          if (H1[i,j,0]==H2[k,m,0] and H1[i,j,1]==H2[k,m,1]):
            for n in range(0,nz):
              if (Om2[k,m,n] > 0 and Om1[i,j,n] > 0):
                Oshift += Om2[k,m,n] - Om1[i,j,n]
                ncommon += 1.0

  Oshift *= (1.0 / ncommon)
  
  # Perform the shifting
  for i in range(0,nx2):
    for j in range(0,ny):
      for k in range(0,nz):
        if (Om2[i,j,k] > 0.0):
          Om2[i,j,k] -= Oshift

  IRmin = nmp.min(H1[:,:,0])
  IRmax = nmp.max(H2[:,:,0])
  dIR = H1[1,0,0] - H1[0,0,0]

  nIRnew = int( (IRmax-IRmin)/dIR ) + 1

  Omnew = nmp.zeros([nIRnew,ny,nz],'d')
  Hnew = nmp.zeros([nIRnew, ny, 2],'d')
  for i in range(0,nIRnew):
    Hnew[i,:,0] = IRmin + dIR*i
    Hnew[i,:,1] = H1[nx1-1,:,1]

  for i in range(0,nIRnew):
    
    ncom = 0
    for j in range(0,nx1):
      if (Hnew[i,0,0] == H1[j,0,0]):
        Omnew[i,:,:] += Om1[j,:,:]
        ncom += 1
        break

    for j in range(0,nx2):
      if (Hnew[i,0,0] == H2[j,0,0]):
        Omnew[i,:,:] += Om2[j,:,:]
        ncom += 1
        break

    if (ncom > 1):
      Omnew[i,:,:] *= 0.5

  return Omnew, Hnew

def write_3D_dos(prefix, Om, H):
  nx = nmp.shape(Om)[0]
  ny = nmp.shape(Om)[1]
  nz = nmp.shape(Om)[2]

  for i in range(0,nz):
    nm = "%s_%d.dat" % (prefix, i)
    otp = open(nm,"w")

    for j in range(0,nx):
    
      for k in range(0,ny):
        line = "%lf %lf %1.10e\n" % (H[j,k,0], H[j,k,1], Om[j,k,i])
        otp.write(line)

      otp.write("\n")

    otp.close()

def get_modA_dimensions(win_id):
  name = "proc%d.input" % (win_id)
  inp = open(name , "r") 
  
  line = inp.readline().split()
  n1 = int( line[2] )

  line = inp.readline().split()
  n2 = int( line[2] )
  
  line = inp.readline().split()
  n3 = int( line[2] )

  inp.close()

  return [n1,n2,n3]

def read_modA_3D_AVG(file_name, n1, n2, n3):
  H = nmp.zeros([n1,n2,n3,3],'d')
  Avg = nmp.zeros([n1,n2,n3],'cdouble')
  Fl = nmp.zeros([n1,n2,n3],'i')


  inp = open( file_name , "r" )
  
  avg = 0.0
  ct = 0
  line_no = 0
  for i in range( 0 , n1 ):
    for j in range( 0 , n2 ):
      for k in range( 0 , n3 ):

        line = inp.readline().split()
        line_no += 1
        while ( len( line ) == 0):
          print i,j,k, line_no
          line = inp.readline().split()
          line_no += 1

        for m in range(0,3):
          H[i,j,k,m] = float( line[m] )
        

        if ( len( line ) == 4 ):
          Fl[i,j,k] = 0
        else:
          Avg[i,j,k] = float( line[3] ) + float( line[4] ) * 1j
          Fl[i,j,k] = 1
          ct += 1

      line = inp.readline()
    line = inp.readline()


  inp.close()

  return [Avg , H, Fl]


def read_modA_3D_DOS(file_name, n1, n2, n3):
  H = nmp.zeros([n1,n2,n3,3],'d')
  Om = nmp.zeros([n1,n2,n3],'d')
  Fl = nmp.zeros([n1,n2,n3],'i')


  inp = open( file_name , "r" )
  
  avg = 0.0
  ct = 0

  for i in range( 0 , n1 ):
    for j in range( 0 , n2 ):
      for k in range( 0 , n3 ):

        line = inp.readline().split()
        for m in range(0,3):
          H[i,j,k,m] = float( line[m] )
        
        Om[i,j,k] = float( line[3] )

        if ( Om[i,j,k] == 0.0 ):
          Fl[i,j,k] = 0
        else:
          Fl[i,j,k] = 1
          avg += Om[i,j,k]
          ct += 1


  inp.close()
  
  return [Om , H, Fl]

# Returns 3D-dos, nQRbins, nIRbins
def read_3D_dos(dir_name, nQIbins):
  
  name = "%s/dos_%d.dat" % (dir_name, 0)
  nIRbins, nQRbins = get_dimensions(name)
 
  Om = nmp.zeros([nIRbins, nQRbins, nQIbins],'d')
  H = nmp.zeros([nIRbins, nQRbins, 2],'d')

  for i in range(0,nQIbins):
    name = "%s/dos_%d.dat" % (dir_name, i)
    Otmp, Htmp = read_dos(name, nIRbins, nQRbins)
    Om[:,:,i] = Otmp
    H[:,:,:] = Htmp

  return [Om, H, nIRbins, nQRbins]




def read_dos(name, nxbins, nybins):
  inp = open(name,"r")
  
  Om = nmp.zeros([nxbins,nybins],'d')
  H = nmp.zeros([nxbins,nybins,2],'d')

  for i in range(0,nxbins):
    for j in range(0,nybins):
      line = inp.readline().split()
      H[i,j,0] = float(line[0])
      H[i,j,1] = float(line[1])
      Om[i,j] = float(line[2])
    inp.readline()

  return Om, H



def get_dimensions(name):
  inp = open(name,"r")

  nQRbins = 0
  nIRbins = 0
 
  line_no = 0
  line = inp.readline().split()
  while (len(line)==3):
    
    line_no += 1
    line = inp.readline().split()
    
    if (len(line)!=3):
      if (nQRbins==0):
        nQRbins = line_no 
      line = inp.readline().split()    
      nIRbins += 1

  inp.close()

  return nIRbins, nQRbins
