#!/usr/bin/env python
import sys
import os
import subprocess
import string
import shutil

print("Test WFEnergyAndForces...")

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
dst1 = 'pseudo.Si'
dst2 = 'pseudo.H'
src1 = sys.argv[nargs-1] + '/' + dst1
src2 = sys.argv[nargs-1] + '/' + dst2

if not os.path.exists(dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)
if not os.path.exists(dst2):
  print("Create link to %s"%dst2)
  os.symlink(src2, dst2)

#run
command = "{} {} -c {} -i {}".format(mpicmd,exe,inp,coords)
print("Run command: {}".format(command))
output = subprocess.check_output(command,shell=True)
lines=output.split(b'\n')

shutil.rmtree('WF')

#analyse output
energies=[]
for line in lines:
  if line.count(b'%%'):
    print(line)
    words=line.split()
    words=words[5].split(b',')[0]
    energy = words.decode()
  if line.count(b'achieved'):
    energies.append(energy)

flag=0
forces=[]
for line in lines:
  if flag>0:
    print(line)
    words=line.split(b'    ')
    forces.append(words[1].decode())
    forces.append(words[2].decode())
    forces.append(words[3].decode())
    flag=flag-1
  if line.count(b'Forces:'):
    flag=2


print("Check energies...")
print( energies )
if len(energies)<2:
  print("Expected two converged energies")
  sys.exit(1)

tol = 1.e-6
diff=eval(energies[1])-eval(energies[0])
print(diff)
if abs(diff)>tol:
  print("Energies differ: {} vs {} !!!".format(energies[0],energies[1]))
  sys.exit(1)

print("Check forces...")
print(forces)
flag=0
for i in range(6):
  diff=eval(forces[i+6])-eval(forces[i])
  print(diff)
  if abs(diff)>1.e-3:
    print("Forces difference larger than tol")
    flag=1
if flag>0:
  sys.exit(1)

print("Test SUCCESSFUL!")
sys.exit(0)
