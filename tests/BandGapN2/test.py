#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test band gap N2...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-5):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-5]
inp = sys.argv[nargs-4]
coords = sys.argv[nargs-3]
print("coordinates file: %s"%coords)

lrs = sys.argv[nargs-2]

#create links to potentials files
dst = 'pseudo.N_ONCVPSP_LDA'
src = sys.argv[-1] + '/' + dst

if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run mgmol
command = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp,coords,lrs)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
#make sure sum forces is below tolerance
lines=output.split(b'\n')

print("Check forces...")
tol = 1.e-5
count = 0
forces = []
for line in lines:
  num_matches = line.count(b'%%')
  if num_matches:
    print(line)
  num_matches = line.count(b'##')
  if num_matches:
    words=line.split()
    if len(words)==8:
      print(line)
      count = count +1
      forces.append(words)

for i in range(5,8):
  print("Force 1: {}".format(forces[1][i]))
  print("Force 2: {}".format(forces[0][i]))
  df=eval(forces[1][i])+eval(forces[0][i])
  print("Force diff. = {}".format(df))
  if abs(df)>tol:
    print("TEST FAILED: Forces not converged!")
    sys.exit(1)

print("Check band gap...")
tole=1.e-1
ExpectedGapEV=8.3
for i in range(len(lines)):
  num_matches1 = lines[i].count(b'convergence')
  num_matches2 = lines[i].count(b'achieved')
  if num_matches1 and num_matches2:
    ii=i+3
    eigenvalues=lines[ii].split()
    print("Eigenvalues: {}".format(eigenvalues) )
    gapHa = eval(eigenvalues[5])-eval(eigenvalues[4])
    gapEv = gapHa*27.211
    if abs(gapEv-ExpectedGapEV)>tole:
      print("ERROR Band gap = {}".format(gapEv))
      print("Expected value: {}".format(ExpectedGapEV))
      sys.exit(1)

sys.exit(0)
