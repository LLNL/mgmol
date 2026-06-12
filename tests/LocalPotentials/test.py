#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test run with local potentials...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-4):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd)) 

exe = sys.argv[nargs-4]
inp = sys.argv[nargs-3]
coords = sys.argv[nargs-2]
print("coordinates file: %s"%coords)

#create links to potentials files
dst = 'pseudo.Li_GTH_PBE'
src = sys.argv[nargs-1] + '/' + dst

cwd = os.getcwd()
if not os.path.exists(cwd+'/'+dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))

output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
#make sure force is below tolerance
lines=output.split(b'\n')

convergence=0
for line in lines:
  if line.count(b'DavidsonSolver') and line.count(b'convergence'):
    convergence=1
    break

if convergence==0:
  print("DavidsonSolver did not converge")
  sys.exit(1)

ended=0
for line in lines:
  if line.count(b'Run') and line.count(b'ended'):
    ended=1

if ended==0:
  print("Run did not end...")
  sys.exit(1)

print("Test PASSED")
sys.exit(0)
