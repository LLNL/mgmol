#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test MVP solver with Chebyshev Approximation ...")

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
dst = 'pseudo.Li_ONCVPSP_LDA'
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

niterations = len(energies)
print("Chebyshev-MVP solver ran for {} iterations".format(niterations))
if niterations>120:
  print("Chebyshev-MVP test FAILED for taking too many iterations")
  sys.exit(1)

energy_ref = -108.264614049136
tol = 1.e-6
print("Check energy...")
last_energy = eval(energies[-1])
print("Energy = {}".format(last_energy))
if abs(last_energy-energy_ref)>tol:
  print("Last energy = {}".format(last_energy))
  sys.exit(1)

tol = 1.e-5
entropy_ref = 0.17653976
print("Check entropy...")
for line in lines:
  if line.count(b'-TS'):
    words=line.split()
    entropy = eval(words[3])
    print("Entropy = {}".format(entropy))
    entropy_diff = entropy+entropy_ref
    if abs(entropy_diff)>tol:
        print("Check entropy test FAILED. Entropy difference = {}".format(abs(entropy_diff)))
        sys.exit(1)
    break

print("Test PASSED")
sys.exit(0)
