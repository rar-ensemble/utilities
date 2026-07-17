import os

def kill_job(name):
  inp = open(name,'r')
  line = inp.readline()
  L = line.split()
  jobID = float(L[5])

  cmd = 'condor_rm %d' % jobID
  os.system(cmd)

  cmd = 'rm %s' % name
  os.system(cmd)


