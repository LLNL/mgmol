#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test Davidson solver with mixing...")

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
dst1 = 'pseudo.H_ONCV_PBE_SG15'
dst2 = 'pseudo.Si_ONCV_PBE_SG15'
src1 = sys.argv[nargs-1] + '/' + dst1
src2 = sys.argv[nargs-1] + '/' + dst2

cwd = os.getcwd()
if not os.path.exists(cwd+'/'+dst1):
  print("Create link to %s"%dst1)
  os.symlink(src1, dst1)
if not os.path.exists(cwd+'/'+dst2):
  print("Create link to %s"%dst2)
  os.symlink(src2, dst2)

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

tol = 2.e-4
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
          print("Error on force component {}".format(i-5))
          sys.exit(1)

niterations = len(energies)
print("Davidson solver ran for {} iterations".format(niterations))
if niterations>50:
  print("Davidson test FAILED for taking too many iterations")
  sys.exit(1)

print("Check energy...")
ref_energy=-6.26479
last_energy = eval(energies[-1])
print("Energy = {}".format(last_energy))
if last_energy>ref_energy:
  print("Last energy = {}".format(last_energy))
  sys.exit(1)

print("Test PASSED")
sys.exit(0)
