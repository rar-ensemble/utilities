import time
import math as m
import numpy as nmp

def sum_Kfile(nm):
  otp = open(nm,'r')
  line = otp.readline().split() 
  sum = 0.0
  while ( len( line ) == 1 ):
    sum += float( line[0] )
    line = otp.readline().split()

  return sum


def generate_seed():
  ch = str( int( time.time() * 10000 ) ) 
  ln = len(ch)
  seed = int( ch[ln-8:ln] )

  return seed

def velocity_line(T):
  seed = generate_seed()
  line = 'velocity all create %f %d mom yes dist gaussian\n' % (T , seed)
  return line
