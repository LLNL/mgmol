#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test Fluorine atom...")

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
dst = 'pseudo.F_ONCV_PBE_SG15'
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

print("Check force...")
tol = 1.e-2
for line in lines:
  num_matches = line.count(b'%%')
  if num_matches:
    print(line)
  num_matches = line.count(b'##')
  if num_matches:
    words=line.split()
    if len(words)==8:
      print(line)
      for i in range(5,8):
        if abs(eval(words[i]))>tol:
          sys.exit(1)

print("Check eigenvalues...")
tole=1.e-3
for i in range(len(lines)):
  num_matches1 = lines[i].count(b'convergence')
  num_matches2 = lines[i].count(b'achieved')
  if num_matches1 and num_matches2:
    ii=i+3
    eigenvalues=lines[ii].split()
    print("Eigenvalues: {}".format(eigenvalues) )
    if abs(eval(eigenvalues[0])+1.070)>tole:
      print("ERROR Eigenvalue 0 = {}".format(eval(eigenvalues[0])))
      sys.exit(1)
    for ii in range(3):
      if abs(eval(eigenvalues[1+ii])+0.410)>tole:
        print("ERROR Eigenvalue {} = {}".format(1+ii,eval(eigenvalues[1+ii])))
        sys.exit(1)
    sys.exit(0)

