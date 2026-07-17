import os
import time
import numpy as nmp
import numpy.linalg as lalg
import random as rand
import math

def explore_min(max_dir):
  # Create and set up the working directories #
  if ( (os.path.exists('find') != 1) or (os.path.exists('2') != 1) ):
    os.system('mkdir find')
    os.system('cp find.input find/simul.input')
    os.system('cp simul.par* find/.')
    
    for i in range(2,max_dir):
      cmd = "mkdir %d" % i
      os.system(cmd)
    
      cmd = "cp simul.par* %d/." % i
      os.system(cmd)
      cmd = "cp follow.input %d/simul.input" % i
      os.system(cmd)
  
  
  # Submits the job that will find the normal modes
  os.system('cp simul.box0 find/.')
  os.system('condor_submit submit-find')
  

  
  # Submits a job to follow the lowest-frequency mode
  os.system('cp submit-blank submit-follow')
  cmd = "echo Initialdir = %d >>submit-follow" % 2
  os.system(cmd)
  cmd = "echo queue >>submit-follow"
  os.system(cmd)
  
  os.system('cp simul.box0 2/.')
  os.system('condor_submit submit-follow | grep cluster >jobID')
  os.system('mv jobID 2/. ;')



  
  # This loop checks on the progress of the mode search.
  # As new modes are found, it submits jobs to follow those
  # modes to saddle points.
  dir = 3
  while (dir < max_dir):
    
    name = './find/mode%d.out' % (dir)
    
    # Wait for the next mode to be printed
    if (os.path.isfile(name) == 0):
      time.sleep(2)
    
    # If a new mode has been found since the last check
    elif (os.path.isfile(name)):

      # Copy all of the mode files into ./dir
      for i in range(3,dir+1):
        cmd = "cp ./find/mode%d.out ./%d/." % (i,dir)
        os.system(cmd)

      cmd = "cp ./simul.box0 ./%d/." % dir
      os.system(cmd)

      os.system('cp submit-blank submit-follow')
      cmd = "echo Initialdir = %d >>submit-follow" % dir
      os.system(cmd)
      cmd = "echo queue >>submit-follow"
      os.system(cmd)
      
      os.system('condor_submit submit-follow | grep cluster >jobID')
      cmd = "mv jobID ./%d/." % dir
      os.system(cmd)

      dir += 1
      
    else:
      print 'What am I doing here??'

  
  # Finished exploring this min!
  print "Victory is mine!"
