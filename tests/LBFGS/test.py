#!/usr/bin/env python
import sys
import os
import subprocess
import string

print("Test LBFGS...")

nargs=len(sys.argv)

mpicmd = sys.argv[1]+" "+sys.argv[2]+" "+sys.argv[3]
for i in range(4,nargs-5):
  mpicmd = mpicmd + " "+sys.argv[i]
print("MPI run command: {}".format(mpicmd))

exe = sys.argv[nargs-5]
inp = sys.argv[nargs-4]
coords = sys.argv[nargs-3]
print("coordinates file: %s"%coords)
lrs = sys.argv[-2]

#create links to potentials files
dst = 'pseudo.Si_ONCV_PBE_SG15'
src = sys.argv[-1] + '/' + dst

if not os.path.exists(dst):
  print("Create link to %s"%dst)
  os.symlink(src, dst)

#run MGmol
command = "{} {} -c {} -i {} -l {}".format(mpicmd,exe,inp,coords,lrs)
output = subprocess.check_output(command,shell=True)

#analyse mgmol standard output
lines=output.split(b'\n')

print("Check LBFGS convergence...")
tol = 4.e-4
force = 1.e6
for line in lines:
  num_matches1 = line.count(b'%%')
  num_matches2 = line.count(b'ION')
  if num_matches1 and num_matches2:
    print(line)

  num_matches1 = line.count(b'max')
  num_matches2 = line.count(b'movable')
  if num_matches1 and num_matches2:
    print(line)
    words=line.split()
    
    force=eval(words[6])

print("Last max force = {}".format(force))

if force>tol:
  sys.exit(1)

sys.exit(0)
