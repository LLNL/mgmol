#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test MVP solver...")

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
dst = 'pseudo.Al_LDA_FHI'
src = sys.argv[nargs-1] + '/' + dst

cwd = os.getcwd()
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

tol = 1.e-4
energies=[]
print("Check forces are smaller than tol = {}".format(tol))
for line in lines:
  if line.count(b'MVP'):
    print(line)
  if line.count(b'%%'):
    print(line)
    words=line.split()
    energy=(words[5].split(b','))[0]
    energies.append(energy)
  if line.count(b'##'):
    words=line.split()
    if len(words)==8:
      print(line)
      for i in range(5,8):
        if abs(eval(words[i]))>tol:
          sys.exit(1)

for line in lines:
  if line.count(b'HDF5-DIAG') and line.count(b'Error'):
    print(line)
    print("Found HDF5 error")
    sys.exit(1)

flag = 0
eigenvalues=[]
for line in lines:
  if line.count(b'FERMI'):
    flag = 0
  if flag==1:
    words=line.split()
    for w in words:
      eigenvalues.append(eval(w))
  if line.count(b'Eigenvalues'):
    flag = 1
    eigenvalues=[]

print(eigenvalues)
tol = 1.e-4
eigenvalue0 = -0.210
if abs(eigenvalues[0]-eigenvalue0)>tol:
  print("Expected eigenvalue 0 to be {}".format(eigenvalue0))
  sys.exit(1)
eigenvalue50 = 0.205
if abs(eigenvalues[50]-eigenvalue50)>tol:
  print("Eeigenvalue 50 = {}".format(eigenvalues[50]))
  print("Expected eigenvalue 50 to be {}".format(eigenvalue50))
  sys.exit(1)

niterations = len(energies)
print("MVP solver ran for {} iterations".format(niterations))
if niterations>180:
  print("MVP test FAILED for taking too many iterations")
  sys.exit(1)

print("Check energy...")
last_energy = eval(energies[-1])
print("Energy = {}".format(last_energy))
if last_energy>-64.390:
  print("Last energy = {}".format(last_energy))
  sys.exit(1)

print("Test PASSED")
sys.exit(0)
