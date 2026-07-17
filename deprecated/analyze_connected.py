import os
import time
import numpy as nmp
import numpy.linalg as lalg
import random as rand
import math as m
import dimw_utils as dimw
import condor

# Waits for searching jobs to finish
def wait(max_dir):
  init_time = time.time()

  # First, must wait for all of the jobs to finish
  # Gives up after two hours have elapsed
  total_time = 0


  # Main waiting loop
  while (condor.shortest_runtime(max_dir) < 3000):
    time.sleep(180)



  # Print condor waiting time info
  final_time = time.time()
  if ( (final_time - init_time) > 600):
    print "I spent %f minutes waiting on Condor" % ( (final_time-init_time)/60)
  else:
    print "I spent %d seconds waiting on Condor!" % ( (final_time-init_time))

  # Kill any remaining jobs
  list = condor.get_running_jobs(max_dir)
  jobs = ' '.join(list)
  if (list!=""):
    cmd = 'condor_rm %s' % jobs
    os.system(cmd)

 # for i in range(2,max_dir):
 #   name1 = '%d/jobID' % i
 #   name2 = '%d/finished.flag' % i
 #   if (os.path.exists(name1) and os.path.exists(name2)==0 ):
 #     condor.kill_job(name1)




Etol = 0.0001
  


# This rountine waits until all saddle searches have finished
# then checks for new connected minima.
# max_dir is the maximum # of modes to be checked
def analyze(cur_state,max_dir, Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE, Temp):
  
  wait(max_dir)
 
  max_new_state = -5

  otp = open('mode_following_traj.dat','a')

  # Loop over directories, read in data file
  for i in range(2, max_dir+1):
    name = '%d/mode_0_trans_data.dat' % i
    
    if (os.path.exists(name)==0):
      continue

    # Reads and tests the validity of the energies
    # and the displacements.  Ensures dimer returned to
    # initial minimum
    flag, swap, E1, E2, Es = dimw.valid_min(i)
    if (flag==0):
      continue

    # Determines rmsd for connected min
    r1, r2, r3 = dimw.min_rmsd(i, swap)
 


    # Checks to see if connected min. is a new minimum
    flag = dimw.new_min(E2, Es, r1, r2, r3, Nstates, DATA, SADDAT)

    if (flag >= 0 and Mcon[flag][cur_state]!=0.0):
      continue
    
    elif (flag>=0 and Mcon[flag][cur_state]==0.0):
      Mcon = dimw.add_connection(cur_state, flag, swap, Mcon, Es, E2, E1, Temp)
      
      line = "%d %d %f %f %f\n" % (cur_state,i,E1,E2,Es)
      otp.write(line)
      
      if (i > max_new_state):
        max_new_state = i
      continue

    else:
      dat = (E1, E2, Es, r1, r2, r3, Temp)
      if (i > max_new_state):
        max_new_state = i
  
      # New, valid minimum has been located.  Add it to our web
      Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE, \
          = dimw.add_state(cur_state,i, swap, Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE, dat)
      
      line = "%d %d %f %f %f\n" % (cur_state,i,E1,E2,Es)
      otp.write(line)

  print "Mode %d returned a new k_{a,b}" % max_new_state
  
  otp.close()

  return Nstates, DATA, STATES, SADDLE, SADDAT, Mcon, NE, 





