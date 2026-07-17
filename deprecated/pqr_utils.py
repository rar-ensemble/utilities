import os, sys
import numpy as np

def write_pqr( name , x , tp , bx , dat1 , dat2 ):
  ns = np.shape(x)[0]


  pqr = open( name , 'w' )
  
  pqr.write('TITLE     Some Junk\n')


  ln = 'CRYST1    %2.3f  %2.3f  %2.3f  90.0  90.0  90.0 P 1        1\n' % \
      (bx[0] , bx[1] , bx[2] )

  pqr.write( ln )
  

  pqr.write('MODEL        1\n')

  fl1 = fl2 = 0

  if ( np.shape(dat1)[0] == ns ):
    fl1 = 1

  if ( np.shape(dat2)[0] == ns ):
    fl2 = 2


  for i in range( 0 , ns ):
    if (fl1):
      if ( fl2 ):
        line = '%s  %d %4s %5s  %4d  %2.5f %2.5f %2.5f %2.4f  %2.4f\n' % \
            ('ATOM' , i , 

  pqr.close()
