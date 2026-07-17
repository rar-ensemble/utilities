import os
import time

def kill_job_list(job_list):
  args = ' '.join(job_list)
  cmd = "condor_rm %s" % args
  os.system(cmd)

def kill_job(name):
  inp = open(name,'r')
  line = inp.readline()
  L = line.split()
  jobID = float(L[5])

  cmd = 'condor_rm %d' % jobID
  os.system(cmd)

  cmd = 'rm %s' % name
  os.system(cmd)

# Checks to see if the condor_q command successfully fetched
# the queue list from the submit node
# Returns 1 if it DID FAIL
def check_failure():
  os.system('grep \"Failed to fetch\" condor-stdout >gp')

  inp = open('gp','r')
  if (inp.readline()==""):
    inp.close()
    return 0
  else:
    inp.close()
    return 1


def time_running(name):
  
  if (os.path.exists(name)==0):
    return 0
    
  inp = open(name,'r')
  line = inp.readline()
  L = line.split()
  jobID = float(L[5])

  
  cmd = 'condor_q %d >condor-stdout' % jobID
  os.system(cmd)
  while (check_failure()):
    time.sleep(300)
    cmd = 'condor_q %d >condor-stdout' % jobID
    os.system(cmd)

  os.system('grep rariggle condor-stdout >py-stdout')

  inp = open('py-stdout','r')
  line = inp.readline()
  if (line==""):
    return 1.0
  L = line.split()
  L = L[4].split(':')
  secs = float(L[2])
  mins = float(L[1])
  L = L[0].split('+')
  hrs = float(L[1])
  
  time = 3600*hrs + 60*mins + secs

  return time

def get_running_jobs(max_dir):
  y = []
  for i in range(2,max_dir):
    name1 = '%d/jobID' % i
    name2 = '%d/finished.flag' % i

    if (os.path.exists(name2)):
      continue

    if (os.path.exists(name1)==0):
      continue
      
    inp = open(name1,'r')
    line = inp.readline()
    L = line.split()
    jobID = float(L[5])

    y.append(str(jobID))

  return y


def shortest_runtime(max_dir):
  min_run_time = 9000.0
  
  y = get_running_jobs(max_dir)

  # If y is NULL, return a small value
  if (y==""):
    return 1.0

  args = ' '.join(y)
  cmd = 'condor_q %s rariggle >condor-stdout' % args
  os.system(cmd)
  
  while (check_failure()):
    time.sleep(300)
    cmd = 'condor_q %s rariggle >condor-stdout' % args
    os.system(cmd)

  os.system('grep rariggle condor-stdout >py-stdout')

  inp = open('py-stdout','r')
  
  kill_list = []
  nkills = 0

  line = "1"
  while (line!=""):
    line = inp.readline()
    if (line==""):
      continue
    L = line.split()
    
    jobID = float(L[0])

    L = L[4].split(':')
    secs = float(L[2])
    mins = float(L[1])
    L = L[0].split('+')
    hrs = float(L[1])
    
    time = 3600*hrs + 60*mins + secs

    if (time < min_run_time):
      min_run_time = time

    if (time > 3600*2.5): #If running for more than 2.5 hours
      kill_list.append(str(jobID))
      nkills += 1

  if (nkills > 0):
    kill_job_list(kill_list)

  inp.close()
  return min_run_time      

