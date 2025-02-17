#!/usr/bin/env python
import sys
import os
import subprocess
import string
import shutil

print("Test DMandEnergyAndForces...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-6):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-5]
inp = sys.argv[nargs-4]
coords = sys.argv[nargs-3]
print("coordinates file: %s"%coords)
lrs = sys.argv[-2]

#create links to potentials files
dst = 'pseudo.N_ONCVPSP_LDA'
src = sys.argv[-1] + '/' + dst

if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run
command = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp,coords,lrs)
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
    break

for line in lines:
  if line.count(b'MVP') and line.count(b'iteration'):
    print(line)
  if line.count(b'Eks2'):
    print(line)
    words=line.split()
    word=words[2]
    energy = word.decode()
    energies.append(energy)
    break

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

print("Test SUCCESSFUL!")
sys.exit(0)
