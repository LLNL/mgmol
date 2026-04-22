#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test computation of dipole moment...")

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
cwd = os.getcwd()

dst = 'pseudo.O_ONCV_PBE_SG15'
src = sys.argv[nargs-1] + '/' + dst
if not os.path.exists(cwd+'/'+dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

dst = 'pseudo.H_ONCV_PBE_SG15'
src = sys.argv[nargs-1] + '/' + dst
if not os.path.exists(cwd+'/'+dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run mgmol
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))

output = subprocess.check_output(command,stderr=subprocess.STDOUT,shell=True)

#analyse mgmol standard output
#make sure force is below tolerance
lines=output.split(b'\n')

convergence=0
for line in lines:
  if line.count(b'DFTsolver:') and line.count(b'convergence'):
    convergence=1
    break

if convergence==0:
  print("MVP Solver did not converge")
  sys.exit(1)

for line in lines:
  if line.count(b'Dipole') and line.count(b'Debye'):
    print(line)
    words=line.split()
    dy = eval(words[4][:-1])
    print("Diploe: {}".format(dy))

tol = 1.e-4
ref_value=2.00519108
if abs(dy-ref_value)>tol:
  print("Expected dipole: {}".format(ref_value))
  sys.exit(1)

print("Test PASSED")
sys.exit(0)
